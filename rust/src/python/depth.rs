use std::collections::HashMap;
use std::collections::HashSet;
use std::path::PathBuf;

use crate::bam;
use crate::cli::DepthOptions;
use crate::io;
use crate::python::utils::{extract_to_option_list, extract_to_option_pathbuf, extract_to_usize};
use pyo3::prelude::*;

#[pymethods]
impl DepthOptions {
    #[new]
    fn new(
        list: Option<HashSet<Vec<u8>>>,
        list_file: Option<PathBuf>,
        bam: Option<PathBuf>,
        cram: Option<PathBuf>,
        fasta: Option<PathBuf>,
        bin_size: usize,
        output: Option<PathBuf>,
    ) -> Self {
        DepthOptions {
            list,
            list_file,
            bam,
            cram,
            fasta,
            bin_size,
            output,
        }
    }
}

#[pyfunction]
pub fn depth_with_options(options: &DepthOptions) -> PyResult<usize> {
    let seq_names = match options.list.to_owned() {
        Some(value) => value,
        _ => match options.list_file.to_owned() {
            value => io::get_list(&value),
        },
    };
    let bam = bam::open_bam(&options.bam, &options.cram, &options.fasta, true);
    bam::get_depth(bam, &seq_names, &options);
    Ok(1)
}

fn convert_hashmap_to_options(py: Python<'_>, map: HashMap<String, PyObject>) -> DepthOptions {
    let list = extract_to_option_list(py, &map, "list");
    let list_file = extract_to_option_pathbuf(py, &map, "list_file");
    let bam = extract_to_option_pathbuf(py, &map, "bam");
    let cram = extract_to_option_pathbuf(py, &map, "cram");
    let fasta = extract_to_option_pathbuf(py, &map, "fasta");
    let output = extract_to_option_pathbuf(py, &map, "output");
    let bin_size = extract_to_usize(py, &map, "bin_size");
    DepthOptions {
        list,
        list_file,
        bam,
        cram,
        fasta,
        bin_size,
        output,
    }
}

#[pyfunction(kwds = "**")]
pub fn depth(py: Python<'_>, kwds: Option<HashMap<String, PyObject>>) -> PyResult<()> {
    let options = match kwds {
        Some(map) => convert_hashmap_to_options(py, map),
        None => panic!["No arguments provided"],
    };
    depth_with_options(&options)?;
    Ok(())
}

// #[pyfunction]
// pub fn depth(py: Python<'_>, map: HashMap<String, PyObject>) -> PyResult<()> {

//     let options = &convert_hashmap_to_options(py, map);
//     depth_with_options(options)?;
//     Ok(())
// }