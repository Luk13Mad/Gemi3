// build with cargo rustc --release -- -C link-arg=-undefined -C link-arg=dynamic_lookup
//rename .dylib to .so mv target/release/libCAVI.dylib Gemi3/CAVI.so
mod CAVI;
use ndarray::Array;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use std::collections::HashMap;
use std::iter::IntoIterator;

#[pyfunction]
fn start_inference(
    verbose: bool,
    debug: bool,
    force: bool,
    MAX_ITERATIONS: u32,
    MAE_THRESHOLD: f64,
    ncgenes: Vec<String>,
    ncguides: Vec<String>,
    sample_names: Vec<String>,
    MAE: Vec<f64>,
    unique_guides: Vec<String>,
    unique_genes: Vec<String>,
    guide_to_gene_map: HashMap<String, String>,
    gene_to_guide_map: HashMap<String, Vec<String>>,
    prio_mu_x: f64,
    prio_tau_x: f64,
    prio_mu_xx: f64,
    prio_tau_xx: f64,
    prio_mu_y: f64,
    prio_tau_y: f64,
    prio_mu_yy: f64,
    prio_tau_yy: f64,
    prio_mu_r: f64,
    prio_tau_r: f64,
    prio_mu_rr: f64,
    prio_tau_rr: f64,
    prio_mu_s: f64,
    prio_tau_s: f64,
    prio_mu_ss: f64,
    prio_tau_ss: f64,
    hash_x: HashMap<String, usize>,
    all_x: Vec<f64>,
    all_xx: Vec<f64>,
    all_x_nrow: usize,
    hash_y: HashMap<String, (usize, usize)>,
    all_y: Vec<f64>,
    all_yy: Vec<f64>,
    all_y_nrow: usize,
    all_y_ncol: usize,
    hash_r: HashMap<String, usize>,
    all_r: Vec<f64>,
    all_rr: Vec<f64>,
    all_r_nrow: usize,
    hash_s: HashMap<String, (usize, usize)>,
    all_s: Vec<f64>,
    all_ss: Vec<f64>,
    all_s_nrow: usize,
    all_s_ncol: usize,
    hash_td: HashMap<String, (usize, usize)>,
    all_tau: Vec<f64>,
    all_tau_nrow: usize,
    all_tau_ncol: usize,
    all_alpha: Vec<f64>,
    all_beta: Vec<f64>,
    all_D: Vec<f64>,
) -> PyResult<Vec<Vec<f64>>> {
    let inference_object = CAVI::CAVI::new(
        verbose,
        debug,
        force,
        MAX_ITERATIONS,
        MAE_THRESHOLD,
        ncgenes,
        ncguides,
        sample_names,
        MAE,
        unique_guides,
        unique_genes,
        guide_to_gene_map,
        gene_to_guide_map,
        prio_mu_x,
        prio_tau_x,
        prio_mu_xx,
        prio_tau_xx,
        prio_mu_y,
        prio_tau_y,
        prio_mu_yy,
        prio_tau_yy,
        prio_mu_r,
        prio_tau_r,
        prio_mu_rr,
        prio_tau_rr,
        prio_mu_s,
        prio_tau_s,
        prio_mu_ss,
        prio_tau_ss,
        hash_x,
        all_x,
        all_xx,
        all_x_nrow,
        hash_y,
        all_y,
        all_yy,
        all_y_nrow,
        all_y_ncol,
        hash_r,
        all_r,
        all_rr,
        all_r_nrow,
        hash_s,
        all_s,
        all_ss,
        all_s_nrow,
        all_s_ncol,
        hash_td,
        all_tau,
        all_tau_nrow,
        all_tau_ncol,
        all_alpha,
        all_beta,
        all_D,
    )
    .unwrap();

    inference_object.start_CAVI().unwrap();

    let out = vec![
        inference_object.MAE.read().unwrap().clone(),
        inference_object.all_x.read().unwrap().clone().to_vec(),
        inference_object.all_xx.read().unwrap().clone().to_vec(),
        inference_object
            .all_y
            .read()
            .unwrap()
            .clone()
            .into_raw_vec(),
        inference_object
            .all_yy
            .read()
            .unwrap()
            .clone()
            .into_raw_vec(),
        inference_object.all_r.read().unwrap().clone().to_vec(),
        inference_object.all_rr.read().unwrap().clone().to_vec(),
        inference_object
            .all_s
            .read()
            .unwrap()
            .clone()
            .into_raw_vec(),
        inference_object
            .all_ss
            .read()
            .unwrap()
            .clone()
            .into_raw_vec(),
        inference_object
            .all_tau
            .read()
            .unwrap()
            .clone()
            .into_raw_vec(),
        inference_object
            .alpha
            .read()
            .unwrap()
            .clone()
            .into_raw_vec(),
        inference_object.beta.read().unwrap().clone().into_raw_vec(),
    ];

    return Ok(out);
}

#[pymodule]
fn rsCAVI(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(start_inference))?;
    Ok(())
}
