use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::ffi::CString;

use blart::TreeMap;

use crate::taxonomy::parse::{Name, Node};
use crate::{taxonomy::parse, utils::styled_progress_bar};

use parse::Nodes;

pub fn build_lookup(
    nodes: &Nodes,
    name_classes: &Vec<String>,
    rank_letter: bool,
) -> HashMap<String, Vec<String>> {
    dbg!(&name_classes);
    let ranks = [
        "subspecies",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
    ];
    let higher_ranks = ["family", "order", "class", "phylum", "kingdom"];
    let mut table = HashMap::new();

    let rank_set: HashSet<&str> = HashSet::from_iter(ranks.iter().cloned());
    let higher_rank_set: HashSet<&str> = HashSet::from_iter(higher_ranks.iter().cloned());
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
    dbg!(&table);
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
    let ranks = [
        "subspecies",
        "species",
        "genus",
        "family",
        // "order",
        // "class",
        // "phylum",
    ];
    let mut matched: HashMap<String, String> = HashMap::new();
    let mut unmatched: HashMap<String, Vec<String>> = HashMap::new();
    let higher_ranks = ["family", "order", "class", "phylum", "kingdom"];
    let higher_rank_set: HashSet<&str> = HashSet::from_iter(higher_ranks.iter().cloned());
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

#[derive(Debug)]
pub struct Ids {
    pub tax_id: String,
    pub rank: String,
    pub anc_ids: HashSet<String>,
}

pub fn build_fast_lookup(nodes: &Nodes, name_classes: &Vec<String>) -> TreeMap<CString, Ids> {
    dbg!(&name_classes);
    let ranks = [
        "subspecies",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom",
    ];
    let higher_ranks = ["family", "order", "class", "phylum", "kingdom"];
    let mut id_map: TreeMap<_, _> = TreeMap::new();

    let rank_set: HashSet<&str> = HashSet::from_iter(ranks.iter().cloned());
    let higher_rank_set: HashSet<&str> = HashSet::from_iter(higher_ranks.iter().cloned());
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
                if name.ends_with("lupus") {
                    dbg!(&name);
                }
                let result = id_map.insert(
                    CString::new(name.clone()).unwrap(),
                    Ids {
                        tax_id: tax_id.clone(),
                        rank: node.rank.clone(),
                        anc_ids: anc_ids.clone(),
                    },
                );
            }
        }
    }

    let to_find = [
        "canis lupus",
        "cnis latrans",
        "canidae",
        "felidae",
        "panthera leo",
    ];
    progress_bar.finish();
    for name in to_find {
        match id_map.get(&CString::new(name).unwrap()) {
            Some(ids) => println!("{name}: {}", ids.tax_id),
            None => {
                println!("{name} is missing");
                let fuzzy: Vec<_> = id_map.fuzzy(&CString::new(name).unwrap(), 2).collect();
                dbg!(fuzzy);
            }
        }
    }

    id_map
}
