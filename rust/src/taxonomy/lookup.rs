use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::ffi::CString;

use blart::TreeMap;

use crate::taxonomy::parse::{Name, Node};
use crate::{taxonomy::parse, utils::styled_progress_bar};

use parse::Nodes;

const RANKS: [&str; 8] = [
    "subspecies",
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "kingdom",
];
const HIGHER_RANKS: [&str; 5] = ["family", "order", "class", "phylum", "kingdom"];

pub fn build_lookup(
    nodes: &Nodes,
    name_classes: &Vec<String>,
    rank_letter: bool,
) -> HashMap<String, Vec<String>> {
    let mut table = HashMap::new();

    let rank_set: HashSet<&str> = HashSet::from_iter(RANKS.iter().cloned());
    let higher_rank_set: HashSet<&str> = HashSet::from_iter(HIGHER_RANKS.iter().cloned());
    let node_count = nodes.nodes.len();
    let progress_bar = styled_progress_bar(node_count, "Building lookup hash");

    for (tax_id, node) in nodes.nodes.iter() {
        progress_bar.inc(1);
        if rank_set.contains(node.rank.as_str()) {
            let lineage = nodes.lineage(&"1".to_string(), tax_id);
            let names = node.names_by_class(Some(&name_classes), true).clone();
            for n in lineage.iter().rev() {
                let n_names = n.names_by_class(Some(&name_classes), true);
                for name in names.iter() {
                    for n_name in n_names.iter() {
                        if higher_rank_set.contains(n.rank.as_str()) {
                            let key = match rank_letter {
                                true => format!(
                                    "{}:{}:{}:{}",
                                    node.rank_letter(),
                                    name,
                                    n.rank_letter(),
                                    n_name
                                ),
                                false => format!("{}:{}", name, n_name),
                            };
                            match table.entry(key) {
                                Entry::Vacant(e) => {
                                    e.insert(vec![node.tax_id()]);
                                }
                                Entry::Occupied(mut e) => {
                                    e.get_mut().push(node.tax_id());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    progress_bar.finish();
    table
}

pub fn build_lineage_lookup(nodes: &Nodes, root_id: &String) -> HashMap<String, String> {
    let node_count = nodes.nodes.len();
    let progress_bar = styled_progress_bar(node_count, "Building lookup hash");
    let mut table = HashMap::new();

    for (tax_id, node) in nodes.nodes.iter() {
        progress_bar.inc(1);
        let lineage = nodes.lineage(root_id, tax_id);
        let s: String = lineage
            .iter()
            .map(|node| node.scientific_name())
            .collect::<Vec<String>>()
            .join("; ");
        let lineage_string = format!("{}; {}; ", s, node.scientific_name());
        table.insert(lineage_string, tax_id.clone());
    }
    progress_bar.finish();
    table
}

pub fn lookup_nodes(
    new_nodes: &Nodes,
    nodes: &mut Nodes,
    new_name_classes: &Vec<String>,
    name_classes: &Vec<String>,
    xref_label: Option<String>,
    create_taxa: bool,
) {
    let mut table = build_lookup(&nodes, &name_classes, true);
    let ranks = RANKS[0..4].to_vec();
    let mut matched: HashMap<String, String> = HashMap::new();
    let mut unmatched: HashMap<String, Vec<String>> = HashMap::new();
    let higher_rank_set: HashSet<&str> = HashSet::from_iter(HIGHER_RANKS.iter().cloned());
    let node_count = new_nodes.nodes.len();
    let progress_bar = styled_progress_bar(node_count, "Looking up names");
    let mut hits = vec![];

    // for (tax_id, node) in new_nodes.nodes.iter() {
    for rank in ranks.into_iter().rev() {
        for node in new_nodes.nodes_by_rank(rank) {
            let tax_id = &node.tax_id;
            progress_bar.inc(1);
            let lineage = new_nodes.lineage(&"1".to_string(), tax_id);
            let names = node.names_by_class(Some(name_classes), true);
            let mut match_tax_id = None;
            let mut hanger_tax_id = None;
            for n in lineage.into_iter().rev() {
                if let Some(match_id) = matched.get(&n.tax_id) {
                    if hanger_tax_id.is_none() {
                        hanger_tax_id = Some(match_id.clone());
                    }
                }
                let n_names = n.names_by_class(Some(new_name_classes), true);
                for name in names.iter() {
                    for n_name in n_names.iter() {
                        if higher_rank_set.contains(n.rank.as_str()) {
                            let key = format!(
                                "{}:{}:{}:{}",
                                node.rank_letter(),
                                name,
                                n.rank_letter(),
                                n_name
                            );
                            match table.get(&key) {
                                None => (),
                                Some(value) => {
                                    if value.len() == 1 {
                                        matched.insert(node.tax_id(), value[0].clone());
                                        match_tax_id = Some(value[0].clone());
                                        break;
                                    }
                                }
                            };
                        }
                    }
                    if match_tax_id.is_some() {
                        break;
                    }
                }
            }
            if let Some(ref_tax_id) = match_tax_id {
                hits.push(ref_tax_id.clone());
                // add node.tax_id to names as an xref
                let names = nodes
                    .nodes
                    .get_mut(&ref_tax_id)
                    .unwrap()
                    .names
                    .as_mut()
                    .unwrap();
                let label = match xref_label {
                    Some(ref l) => l.clone(),
                    None => "".to_string(),
                };
                names.push(Name {
                    tax_id: ref_tax_id.clone(),
                    name: node.tax_id(),
                    unique_name: format!("{}:{}", &label, node.tax_id()),
                    class: xref_label.clone(),
                });
                break;
            } else if create_taxa {
                if let Some(hanger_id) = hanger_tax_id {
                    // Create new node and hang on hanger_tax_id
                    let new_tax_id = match xref_label {
                        Some(ref l) => format!("{}:{}", l, node.tax_id()),
                        None => format!(":{}", node.tax_id()),
                    };
                    matched.insert(node.tax_id(), new_tax_id.clone());

                    nodes.nodes.insert(
                        new_tax_id.clone(),
                        Node {
                            tax_id: new_tax_id.clone(),
                            parent_tax_id: hanger_id.clone(),
                            names: match node.names.clone() {
                                Some(names) => Some(
                                    names
                                        .iter()
                                        .map(|n| Name {
                                            tax_id: new_tax_id.clone(),
                                            ..n.clone()
                                        })
                                        .collect(),
                                ),
                                None => None,
                            },
                            rank: node.rank(),
                            scientific_name: node.scientific_name.clone(),
                        },
                    );
                    match nodes.children.entry(hanger_id.clone()) {
                        Entry::Vacant(e) => {
                            e.insert(vec![new_tax_id.clone()]);
                        }
                        Entry::Occupied(mut e) => {
                            e.get_mut().push(new_tax_id.clone());
                        }
                    }
                    let parent_node = nodes.nodes.get(&hanger_id).unwrap();
                    let key = format!(
                        "{}:{}:{}:{}",
                        node.rank_letter(),
                        node.lc_scientific_name(),
                        parent_node.rank_letter(),
                        parent_node.lc_scientific_name()
                    );
                    match table.entry(key) {
                        Entry::Vacant(e) => {
                            e.insert(vec![new_tax_id]);
                        }
                        Entry::Occupied(mut e) => {
                            e.get_mut().push(new_tax_id);
                        }
                    }
                } else {
                    match unmatched.entry(node.rank()) {
                        Entry::Vacant(e) => {
                            e.insert(vec![node.lc_tax_id()]);
                        }
                        Entry::Occupied(mut e) => {
                            e.get_mut().push(node.lc_tax_id());
                        }
                    }
                }
            }
        }
    }
    progress_bar.finish();
    // for rank in ranks {
    //     eprintln!(
    //         "{:?}: {:?}, {:?}",
    //         rank,
    //         match matched.entry(rank.to_string()) {
    //             Entry::Vacant(_) => 0,
    //             Entry::Occupied(e) => {
    //                 e.get().len()
    //             }
    //         },
    //         match unmatched.entry(rank.to_string()) {
    //             Entry::Vacant(_) => 0,
    //             Entry::Occupied(e) => {
    //                 e.get().len()
    //             }
    //         },
    //     )
    // }
}

#[derive(Clone, Debug, Default)]
pub struct TaxonInfo {
    pub tax_id: String,
    pub name: String,
    pub rank: String,
    pub anc_ids: HashSet<String>,
}

pub fn build_fast_lookup(
    nodes: &Nodes,
    name_classes: &Vec<String>,
) -> TreeMap<CString, Vec<TaxonInfo>> {
    let mut id_map: TreeMap<_, _> = TreeMap::new();

    let rank_set: HashSet<&str> = HashSet::from_iter(RANKS.iter().cloned());
    let higher_rank_set: HashSet<&str> = HashSet::from_iter(HIGHER_RANKS.iter().cloned());
    let node_count = nodes.nodes.len();
    let progress_bar = styled_progress_bar(node_count, "Building lookup hash");

    for (tax_id, node) in nodes.nodes.iter() {
        progress_bar.inc(1);
        if rank_set.contains(node.rank.as_str()) {
            let lineage = nodes.lineage(&"1".to_string(), tax_id);
            let names = node.names_by_class(Some(&name_classes), true);
            let anc_ids: HashSet<String> = lineage
                .iter()
                .filter(|n| higher_rank_set.contains(n.rank.as_str()))
                .map(|n| n.tax_id())
                .collect();
            for name in names {
                let key = name.clone();
                let taxon_info = TaxonInfo {
                    tax_id: tax_id.clone(),
                    name: node.scientific_name(),
                    rank: node.rank(),
                    anc_ids: anc_ids.clone(),
                };
                match id_map.entry(CString::new(key.clone()).unwrap()) {
                    blart::map::Entry::Vacant(e) => {
                        e.insert(vec![taxon_info]);
                    }
                    blart::map::Entry::Occupied(mut e) => {
                        e.get_mut().push(taxon_info);
                    }
                }
            }
        }
    }

    progress_bar.finish();

    id_map
}

#[derive(Clone, Debug, Default)]
pub struct Candidate {
    pub name: String,
    pub tax_id: Option<String>,
    pub rank: String,
    pub anc_ids: Option<HashSet<String>>,
}

#[derive(Debug, Clone, Default)]
pub enum MatchStatus {
    Match(Candidate),
    Mismatch(Candidate),
    MultiMatch(Vec<Candidate>),
    PutativeMatch(Candidate),
    #[default]
    None,
}

#[derive(Debug, Clone, Default)]
pub struct TaxonMatch {
    pub taxon: Candidate,
    pub taxon_id: Option<String>,
    pub rank_status: Option<MatchStatus>,
    pub rank_options: Option<Vec<Candidate>>,
    pub higher_status: Option<MatchStatus>,
    pub higher_options: Option<Vec<Candidate>>,
}

fn check_higher_taxon(taxon: &Candidate, higher_taxon: &Candidate) {
    let higher_tax_id = higher_taxon.clone().tax_id.unwrap();
    if taxon.anc_ids.clone().unwrap().contains(&higher_tax_id) {
        println!("Taxon {} has taxID {}", higher_taxon.name, higher_tax_id);
    } else {
        println!(
            "Taxon {} is not an ancestor of {}",
            taxon.name,
            taxon.tax_id.clone().unwrap(),
        );
    }
}

fn check_higher_rank(taxon: &Candidate, taxon_match: &TaxonMatch) {
    match taxon_match.higher_status.clone() {
        Some(MatchStatus::Match(_)) => (),
        Some(MatchStatus::Mismatch(_)) => (),
        Some(MatchStatus::MultiMatch(higher_taxa)) => {
            println!("Taxon {} has multiple matches", taxon_match.taxon.name);
        }
        Some(MatchStatus::PutativeMatch(higher_taxon)) => {
            check_higher_taxon(taxon, &higher_taxon);
        }
        _ => {
            if let Some(higher_options) = taxon_match.higher_options.clone() {
                for higher_taxon in higher_options.iter() {
                    check_higher_taxon(taxon, &higher_taxon);
                }
            }
            // println!("No match for taxon name {}", taxon_match.taxon.name);
        }
    }
}

pub fn match_taxonomy_section(
    taxonomy_section: &HashMap<String, String>,
    id_map: &TreeMap<CString, Vec<TaxonInfo>>,
) -> TaxonMatch {
    let mut taxon_id = taxonomy_section.get("taxon_id");
    if taxon_id.is_none() {
        //
    } else if let Some(tax_id) = taxon_id {
        if tax_id == "None" {
            taxon_id = None;
        }
    } else if let Some(ids) = id_map.get(&CString::new(taxon_id.unwrap().clone()).unwrap()) {
        if ids.len() == 1 {
            return TaxonMatch {
                taxon: Candidate {
                    tax_id: Some(ids[0].tax_id.clone()),
                    rank: ids[0].rank.clone(),
                    name: ids[0].name.clone(),
                    anc_ids: Some(ids[0].anc_ids.clone()),
                },
                taxon_id: Some(ids[0].tax_id.clone()),
                rank_status: Some(MatchStatus::Match(Candidate {
                    tax_id: Some(ids[0].tax_id.clone()),
                    rank: ids[0].rank.clone(),
                    name: ids[0].name.clone(),
                    anc_ids: Some(ids[0].anc_ids.clone()),
                })),
                ..Default::default()
            };
        }
    }

    let mut ranks = vec![];
    if taxonomy_section.contains_key("taxon") {
        ranks.push("taxon".to_string());
    }
    for rank in RANKS.iter() {
        if taxonomy_section.contains_key(*rank) {
            ranks.push(rank.to_string());
        }
    }
    let lower_ranks: HashSet<&str> = RANKS[0..3].iter().cloned().collect();

    let mut taxon_match = TaxonMatch::default();
    for (i, rank) in ranks.iter().enumerate() {
        let mut name = taxonomy_section.get(rank).unwrap().clone();
        let taxon = Candidate {
            name: name.clone(),
            tax_id: taxon_id.clone().map(|s| s.clone()),
            rank: rank.clone(),
            ..Default::default()
        };
        if i == 0 {
            taxon_match = TaxonMatch {
                taxon: taxon.clone(),
                ..Default::default()
            };
        } else if lower_ranks.contains(rank.as_str()) {
            continue;
        }
        name = name
            .chars()
            .map(|c| {
                if c.is_alphanumeric() {
                    c.to_ascii_lowercase()
                } else {
                    ' '
                }
            })
            .collect();
        match id_map.get(&CString::new(name.clone()).unwrap()) {
            Some(ids) => {
                if ids.len() > 1 {
                    let mut candidates = vec![];
                    for id in ids.iter() {
                        candidates.push(Candidate {
                            tax_id: Some(id.tax_id.clone()),
                            rank: id.rank.clone(),
                            name: id.name.clone(),
                            anc_ids: Some(id.anc_ids.clone()),
                        });
                    }
                    if i == 0 {
                        if let Some(tax_id) = taxon_match.clone().taxon.tax_id {
                            let mut has_match = false;
                            for candidate in candidates.iter() {
                                if tax_id == candidate.tax_id.clone().unwrap() {
                                    taxon_match.rank_status =
                                        Some(MatchStatus::Match(taxon.clone()));
                                    taxon_match.taxon_id = Some(candidate.tax_id.clone().unwrap());
                                    has_match = true;
                                    break;
                                }
                            }
                            if !has_match {
                                println!(
                                    "multi mismatch for {}: {} at rank {}",
                                    taxon.name, tax_id, &rank
                                );
                                // if i == 0 {
                                //     taxon_match.rank_status = Some(MatchStatus::MultiMatch(candidates));
                                // } else {
                                //     taxon_match.higher_status =
                                //         Some(MatchStatus::MultiMatch(candidates));
                                // }
                            }
                        } else {
                            taxon_match.rank_status = Some(MatchStatus::MultiMatch(candidates));
                        }
                    } else {
                        taxon_match.higher_status = Some(MatchStatus::MultiMatch(candidates));
                    }
                } else {
                    let ids = ids.first().unwrap();
                    if i == 0 {
                        // Same rank as record
                        if let Some(tax_id) = taxon_id {
                            // has taxon ID
                            if tax_id.clone() == ids.tax_id {
                                // Exact match
                                taxon_match.rank_status = Some(MatchStatus::Match(taxon.clone()));
                                taxon_match.taxon_id = Some(ids.tax_id.clone());
                                break;
                            } else {
                                // Mismatched taxon_id, possible namespace collision
                                taxon_match.rank_status = Some(MatchStatus::Mismatch(Candidate {
                                    tax_id: Some(ids.tax_id.clone()),
                                    ..taxon.clone()
                                }));
                            }
                        } else {
                            // no taxon ID
                            taxon_match.rank_status = Some(MatchStatus::PutativeMatch(Candidate {
                                tax_id: Some(ids.tax_id.clone()),
                                anc_ids: Some(ids.anc_ids.clone()),
                                rank: ids.rank.clone(),
                                name: ids.name.clone(),
                            }));
                        }
                    } else {
                        // Higher rank
                        // if anc_ids.contains(&ids.tax_id) {
                        //     taxon_match = TaxonMatch {
                        //         higher_status: Some(MatchStatus::Match(Candidate {
                        //             tax_id: Some(ids.tax_id.clone()),
                        //             ..taxon.clone()
                        //         })),
                        //         ..taxon_match
                        //     };
                        // } else
                        // if anc_ids.len() > 0 {
                        //     taxon_match = TaxonMatch {
                        //         higher_status: Some(MatchStatus::Mismatch(Candidate {
                        //             tax_id: Some(ids.tax_id.clone()),
                        //             rank: ids.rank.clone(),
                        //             name: ids.name.clone(),
                        //             anc_ids: Some(ids.anc_ids.clone()),
                        //         })),
                        //         ..taxon_match
                        //     };
                        // } else {
                        taxon_match = TaxonMatch {
                            higher_status: Some(MatchStatus::PutativeMatch(Candidate {
                                tax_id: Some(ids.tax_id.clone()),
                                rank: ids.rank.clone(),
                                name: ids.name.clone(),
                                anc_ids: Some(ids.anc_ids.clone()),
                            })),
                            ..taxon_match
                        };
                        // }
                        break;
                    }
                }
            }
            None => {
                // println!("{name} is missing");
                let fuzzy: Vec<_> = id_map
                    .fuzzy(&CString::new(name.clone()).unwrap(), 2)
                    .collect();
                if fuzzy.len() > 0 {
                    let mut candidates = vec![];
                    for fuzzies in fuzzy.iter() {
                        for f in fuzzies.1.iter() {
                            if i > 0 || f.rank == taxon_match.taxon.rank {
                                candidates.push(Candidate {
                                    tax_id: Some(f.tax_id.clone()),
                                    rank: f.rank.clone(),
                                    name: f.name.clone(),
                                    anc_ids: Some(f.anc_ids.clone()),
                                });
                            }
                        }
                        // if i > 0 || f.1.rank == taxon_match.taxon.rank {
                        //     candidates.push(Candidate {
                        //         tax_id: Some(f.1.tax_id.clone()),
                        //         rank: f.1.rank.clone(),
                        //         name: f.1.name.clone(),
                        //         anc_ids: Some(f.1.anc_ids.clone()),
                        //     });
                        // }
                    }
                    if candidates.len() > 0 {
                        if i == 0 {
                            taxon_match.rank_options = Some(candidates);
                        } else {
                            taxon_match.higher_options = Some(candidates);
                        }
                    }
                }
            }
        }
    }
    match taxon_match.rank_status.clone() {
        Some(MatchStatus::Match(taxon)) => {
            // println!("Taxon {} has taxID {}", taxon.name, taxon.tax_id.unwrap());
        }
        Some(MatchStatus::Mismatch(taxon)) => {
            // println!(
            //     "Taxon {} has mismatched taxID, {} != {}",
            //     taxon_match.taxon.name,
            //     taxon_match.taxon.tax_id.clone().unwrap(),
            //     taxon.tax_id.unwrap()
            // );
            // TODO: check merged IDs
        }
        Some(MatchStatus::MultiMatch(taxa)) => {
            println!("Taxon {} has multiple matches", taxon_match.taxon.name);
            for taxon in taxa.iter() {
                println!(
                    "checking match to {}: {}",
                    taxon.name,
                    taxon.clone().tax_id.unwrap()
                );
                check_higher_rank(&taxon, &taxon_match);
            }
        }
        Some(MatchStatus::PutativeMatch(taxon)) => {
            // println!(
            //     "Taxon {} has putative match to {}",
            //     taxon_match.taxon.name,
            //     taxon.clone().tax_id.unwrap()
            // );
            // check_higher_rank(&taxon, &taxon_match);
        }
        _ => {
            if let Some(rank_options) = taxon_match.rank_options.clone() {
                for taxon in rank_options.iter() {
                    // println!(
                    //     "Taxon {} has potential match to {}, {}",
                    //     taxon_match.taxon.name,
                    //     taxon.name,
                    //     taxon.tax_id.clone().unwrap()
                    // );
                    // check_higher_rank(&taxon, &taxon_match);
                }
            }
            // println!("No match for taxon name {}", taxon_match.taxon.name);
        }
    }
    taxon_match
}
