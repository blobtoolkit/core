//!
//! Invoked by calling:
//! `blobtk taxonomy <args>`

use std::collections::HashMap;
use std::collections::HashSet;

use anyhow;
use fst::{IntoStreamer, Set, Streamer};

// use std::time::{Duration, Instant};

use crate::cli;
use crate::error;
use crate::io;

/// Functions for ncbi taxonomy processing.
pub mod parse;

/// Functions for name lookup.
pub mod lookup;

pub use cli::TaxonomyOptions;

pub use lookup::{build_lookup, lookup_nodes};

use self::lookup::build_fuzzy_lookup;
use self::parse::Name;
use self::parse::{
    parse_ena_jsonl, parse_file, parse_gbif, parse_taxdump, write_taxdump, Node, Nodes,
};

// use std::error::Error;
// use csv::Reader;

// fn example() -> Result<(), Box<dyn Error>> {
//     let mut rdr = Reader::from_path("foo.csv")?;
//     for result in rdr.records() {
//         let record = result?;
//         println!("{:?}", record);
//     }
//     Ok(())
// }

fn load_options(options: &cli::TaxonomyOptions) -> Result<cli::TaxonomyOptions, error::Error> {
    if let Some(config_file) = options.config_file.clone() {
        let reader = match io::file_reader(config_file.clone()) {
            Some(r) => r,
            None => {
                return Err(error::Error::FileNotFound(format!(
                    "{}",
                    &config_file.to_str().unwrap()
                )))
            }
        };
        let taxonomy_options: cli::TaxonomyOptions = match serde_yaml::from_reader(reader) {
            Ok(options) => options,
            Err(err) => {
                return Err(error::Error::SerdeError(format!(
                    "{} {}",
                    &config_file.to_str().unwrap(),
                    err.to_string()
                )))
            }
        };
        return Ok(TaxonomyOptions {
            path: match taxonomy_options.path {
                Some(path) => Some(path),
                None => options.path.clone(),
            },
            taxonomy_format: match taxonomy_options.taxonomy_format {
                Some(taxonomy_format) => Some(taxonomy_format),
                None => options.taxonomy_format.clone(),
            },
            root_taxon_id: match taxonomy_options.root_taxon_id {
                Some(root_taxon_id) => Some(root_taxon_id),
                None => options.root_taxon_id.clone(),
            },
            base_taxon_id: match taxonomy_options.base_taxon_id {
                Some(base_taxon_id) => Some(base_taxon_id),
                None => options.base_taxon_id.clone(),
            },
            out: match taxonomy_options.out {
                Some(out) => Some(out),
                None => options.out.clone(),
            },
            xref_label: match taxonomy_options.xref_label {
                Some(xref_label) => Some(xref_label),
                None => options.xref_label.clone(),
            },
            name_classes: if taxonomy_options.name_classes.len() > 0 {
                taxonomy_options.name_classes.clone()
            } else {
                options.name_classes.clone()
            },
            create_taxa: taxonomy_options.create_taxa.clone(),
            taxonomies: taxonomy_options.taxonomies.clone(),
            genomehubs_files: match taxonomy_options.genomehubs_files {
                Some(genomehubs_files) => Some(genomehubs_files),
                None => options.genomehubs_files.clone(),
            },

            ..Default::default()
        });
    }
    Ok(options.clone())
}

fn taxdump_to_nodes(
    options: &cli::TaxonomyOptions,
    existing: Option<&mut Nodes>,
) -> Result<Nodes, error::Error> {
    let options = load_options(&options)?;
    let nodes;
    if let Some(taxdump) = options.path.clone() {
        nodes = match options.taxonomy_format {
            Some(cli::TaxonomyFormat::NCBI) => parse_taxdump(taxdump).unwrap(),
            Some(cli::TaxonomyFormat::GBIF) => parse_gbif(taxdump).unwrap(),
            Some(cli::TaxonomyFormat::ENA) => parse_ena_jsonl(taxdump, existing).unwrap(),
            None => {
                return Err(error::Error::FileNotFound(format!(
                    "{}",
                    &taxdump.to_str().unwrap()
                )))
            }
        };
    } else {
        return Err(error::Error::NotDefined(format!("taxdump")));
    }
    Ok(nodes)
}

fn get_ranks_from_row(row: HashMap<String, HashMap<String, String>>) -> Vec<String> {
    let mut ranks = vec![];
    let wanted_ranks = vec![
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ];

    if row.len() > 0 {
        let row_taxonomy = row.get("taxonomy").unwrap();
        for rank in &wanted_ranks {
            if row_taxonomy.get(*rank).is_some() {
                ranks.push(rank.to_string());
            }
        }
    }

    ranks
}

fn extract_ranks(
    taxonomy: &HashMap<String, String>,
    ranks: &Vec<String>,
) -> (HashMap<String, String>, String) {
    let mut extracted_ranks = HashMap::new();
    let mut lowest = "".to_string();
    for rank in ranks {
        if let Some(rank_value) = taxonomy.get(rank) {
            extracted_ranks.insert(rank.to_string(), rank_value.clone());
            lowest = rank_value.clone();
        }
    }
    (extracted_ranks, lowest)
}

fn lookup_rows(
    nodes: &mut Nodes,
    rows: &mut Vec<HashMap<String, HashMap<String, String>>>,
    table: &mut HashMap<String, Vec<String>>,
    fuzzy_table: Option<&Set<Vec<u8>>>,
) -> (
    Nodes,
    Vec<HashMap<String, HashMap<String, String>>>,
    Vec<HashMap<String, HashMap<String, String>>>,
    Vec<HashMap<String, Vec<String>>>,
) {
    let ranks = get_ranks_from_row(rows[0].clone());

    let mut matched_nodes = HashMap::new();
    let mut matched_ids = HashSet::new();
    let mut novel_ids = HashSet::new();
    let mut matched_children = HashMap::new();
    let mut matched_rows = vec![];
    let mut unmatched_rows = vec![];
    let mut spellings = vec![];
    for row in rows.iter_mut() {
        let mut row_taxonomy = row.get_mut("taxonomy").unwrap();
        let taxon_id = row_taxonomy.get("taxon_id");
        if taxon_id.is_some()
            && taxon_id != Some(&"None".to_string())
            && taxon_id != Some(&"NA".to_string())
            && taxon_id != Some(&"".to_string())
        {
            if let Some(node) = nodes.nodes.get(taxon_id.unwrap()) {
                matched_rows.push(row.clone());
                matched_nodes.insert(node.tax_id(), node.clone());
                if let Some(children) = nodes.children.get(node.tax_id.as_str()) {
                    matched_children.insert(node.tax_id(), children.clone());
                }
                let lineage = nodes.lineage(&"root".to_string(), &node.tax_id());
                for anc_node in lineage {
                    if let None = matched_ids.get(anc_node.tax_id.as_str()) {
                        matched_ids.insert(anc_node.tax_id.clone());
                    }
                }
            } else {
                unmatched_rows.push(row.clone());
            }
        } else {
            let (extracted_ranks, lowest_rank) = extract_ranks(row_taxonomy, &ranks);
            let mut previous_rank_value = "root".to_string();

            let mut row_nodes = HashMap::new();
            let row_children = HashMap::new();
            for rank in &ranks {
                if let Some(rank_value) = extracted_ranks.get(rank) {
                    let names = vec![Name {
                        name: rank_value.clone(),
                        class: Some("scientific name".to_string()),
                        ..Default::default()
                    }];
                    let new_node = Node {
                        tax_id: rank_value.clone(),
                        parent_tax_id: previous_rank_value.clone(),
                        rank: rank.to_string(),
                        scientific_name: Some(rank_value.clone()),
                        names: Some(names),
                        ..Default::default()
                    };
                    row_nodes.insert(rank_value.clone(), new_node.clone());
                    novel_ids.insert(rank_value.clone());
                    previous_rank_value = rank_value.clone();
                }
            }

            let new_nodes = Nodes {
                nodes: row_nodes,
                children: row_children,
            };
            let (matched, spellcheck) = lookup_nodes(
                &new_nodes,
                nodes,
                table,
                fuzzy_table,
                &vec!["scientific name".to_string()],
                &vec!["scientific name".to_string()],
                None,
                false,
            );
            if let Some(taxid) = matched.get(&lowest_rank) {
                row_taxonomy.insert("taxon_id".to_string(), taxid.clone());
                matched_rows.push(row.clone());
                if let None = matched_ids.get(taxid) {
                    matched_ids.insert(taxid.clone());
                    matched_nodes.insert(taxid.clone(), nodes.nodes.get(taxid).unwrap().clone());
                    if let Some(children) = nodes.children.get(taxid) {
                        matched_children.insert(taxid.clone(), children.clone());
                    }
                }
            } else {
                unmatched_rows.push(row.clone());
                if spellcheck.len() > 0 {
                    spellings.push(spellcheck);
                }
            }
            // dbg!(new_nodes);

            // pub struct Node {
            //     pub tax_id: String,
            //     pub parent_tax_id: String,
            //     pub rank: String,
            //     pub names: Option<Vec<Name>>,
            //     pub scientific_name: Option<String>,
            // }
        }
    }

    (
        Nodes {
            nodes: matched_nodes,
            children: matched_children,
        },
        matched_rows,
        unmatched_rows,
        spellings,
    )
}

/// Execute the `taxonomy` subcommand from `blobtk`.
pub fn taxonomy(options: &cli::TaxonomyOptions) -> Result<(), anyhow::Error> {
    let options = load_options(&options)?;
    let mut nodes = taxdump_to_nodes(&options, None).unwrap();
    // if let Some(taxdump) = options.path.clone() {
    //     nodes = match options.taxonomy_format {
    //         Some(cli::TaxonomyFormat::NCBI) => parse_taxdump(taxdump)?,
    //         Some(cli::TaxonomyFormat::GBIF) => parse_gbif(taxdump)?,
    //         None => Nodes {
    //             ..Default::default()
    //         },
    //     };
    //     if let Some(taxdump_out) = options.out.clone() {
    //         let root_taxon_ids = options.root_taxon_id.clone();
    //         let base_taxon_id = options.base_taxon_id.clone();
    //         write_taxdump(&nodes, root_taxon_ids, base_taxon_id, taxdump_out);
    //     }
    // }

    if let Some(taxonomies) = options.taxonomies.clone() {
        for taxonomy in taxonomies {
            let new_nodes = taxdump_to_nodes(&taxonomy, Some(&mut nodes)).unwrap();
            // match new_nodes to nodes
            if let Some(taxonomy_format) = taxonomy.taxonomy_format {
                if matches!(taxonomy_format, cli::TaxonomyFormat::ENA) {
                    continue;
                }
                let mut table = build_lookup(&nodes, &options.name_classes, true);
                lookup_nodes(
                    &new_nodes,
                    &mut nodes,
                    &mut table,
                    None,
                    &taxonomy.name_classes,
                    &options.name_classes,
                    taxonomy.xref_label.clone(),
                    taxonomy.create_taxa,
                );
            }
        }
    }

    if let Some(genomehubs_files) = options.genomehubs_files.clone() {
        // dbg!(&options);
        let mut table = build_lookup(&nodes, &options.name_classes, true);
        let fuzzy_table = build_fuzzy_lookup(&nodes, &options.name_classes, true);
        let mut all_nodes = Nodes {
            nodes: HashMap::new(),
            children: HashMap::new(),
        };
        for genomehubs_file in genomehubs_files {
            // match taxa to nodes
            let mut rows = parse_file(genomehubs_file, &table)?;
            let (matched_nodes, matched_rows, unmatched_rows, spellings) =
                lookup_rows(&mut nodes, &mut rows, &mut table, Some(&fuzzy_table));
            dbg!(rows.len());
            dbg!(matched_rows.len());
            dbg!(unmatched_rows.len());
            all_nodes.merge(&matched_nodes, &nodes);
        }

        if let Some(taxdump_out) = options.out.clone() {
            let root_taxon_ids = options.root_taxon_id.clone();
            let base_taxon_id = options.base_taxon_id.clone();
            write_taxdump(&all_nodes, root_taxon_ids, base_taxon_id, taxdump_out);
        }
    }

    // if let Some(taxdump_out) = options.out.clone() {
    //     let root_taxon_ids = options.root_taxon_id.clone();
    //     let base_taxon_id = options.base_taxon_id.clone();
    //     write_taxdump(&nodes, root_taxon_ids, base_taxon_id, taxdump_out);
    // }

    // if let Some(gbif_backbone) = options.gbif_backbone.clone() {
    //     // let trie = build_trie(&nodes);
    //     if let Ok(gbif_nodes) = parse_gbif(gbif_backbone) {
    //         println!("{}", gbif_nodes.nodes.len());
    //         if let Some(taxdump_out) = options.taxdump_out.clone() {
    //             let root_taxon_ids = options.root_taxon_id.clone();
    //             let base_taxon_id = options.base_taxon_id.clone();
    //             write_taxdump(&gbif_nodes, root_taxon_ids, base_taxon_id, taxdump_out);
    //         }
    //     }
    // }

    // if let Some(data_dir) = options.data_dir.clone() {
    //     let trie = build_trie(&nodes);
    //     let rank = "genus".to_string();
    //     let higher_rank = "family".to_string();
    //     let start = Instant::now();
    //     dbg!(trie.predictive_search(vec![
    //         rank,
    //         "arabidopsis".to_string(),
    //         higher_rank,
    //         "brassicaceae".to_string()
    //     ]));
    //     let duration = start.elapsed();

    //     println!("Time elapsed in expensive_function() is: {:?}", duration);
    // }
    // TODO: make lookup case insensitive
    // TODO: add support for synonym matching
    // TODO: read in taxon names from additonal files
    // TODO: add support for fuzzy matching?
    // TODO: hang additional taxa on the loaded taxonomy
    Ok(())
}
