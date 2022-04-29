#![allow(non_snake_case)]
#![allow(non_camel_case_types)]
#![allow(dead_code)]
use itertools::Itertools;
use itertools::{iproduct, izip};
use ndarray::ShapeBuilder;
use ndarray::{Array, Ix1, Ix2};
use std::collections::HashMap;
use std::sync::RwLock;
mod update_x;
mod update_y;
mod update_r;
mod update_s;
mod update_tau;
mod update_MAE;
mod lookup_functions;
mod utility;

#[derive(Debug)]
pub struct CAVI<'a> {
    pub verbose: bool,
    pub debug: bool,
    pub force: bool,
    pub MAX_ITERATIONS: u32,
    pub MAE_THRESHOLD: f64,
    pub prio_mu_x: f64, //x and xx prio
    pub prio_tau_x: f64,
    pub prio_mu_xx: f64,
    pub prio_tau_xx: f64,
    pub prio_mu_y: f64, // y and yy prio
    pub prio_tau_y: f64,
    pub prio_mu_yy: f64,
    pub prio_tau_yy: f64,
    pub prio_mu_r: f64, // r and rr prio
    pub prio_tau_r: f64,
    pub prio_mu_rr: f64,
    pub prio_tau_rr: f64,
    pub prio_mu_s: f64, // s and ss prio
    pub prio_tau_s: f64,
    pub prio_mu_ss: f64,
    pub prio_tau_ss: f64,
    pub ncgenes: Vec<String>,      //ncgenes
    pub ncguides: Vec<String>,     //ncguides. Guides corresponding to ncgenes
    pub sample_names: Vec<String>, //sample names
    pub MAE: RwLock<Vec<f64>>,
    pub hash_x: HashMap<String, usize>,
    pub all_x: RwLock<ndarray::Array<f64, ndarray::Ix1>>,
    pub all_xx: RwLock<ndarray::Array<f64, ndarray::Ix1>>,
    pub hash_y: HashMap<String, (usize, usize)>,
    pub all_y: RwLock<ndarray::Array<f64, ndarray::Ix2>>,
    pub all_yy: RwLock<ndarray::Array<f64, ndarray::Ix2>>,
    pub hash_r: HashMap<String, usize>,
    pub all_r: RwLock<ndarray::Array<f64, ndarray::Ix1>>,
    pub all_rr: RwLock<ndarray::Array<f64, ndarray::Ix1>>,
    pub hash_s: HashMap<String, (usize, usize)>,
    pub all_s: RwLock<ndarray::Array<f64, ndarray::Ix2>>,
    pub all_ss: RwLock<ndarray::Array<f64, ndarray::Ix2>>,
    pub hash_td: HashMap<String, (usize, usize)>,
    pub all_tau: RwLock<ndarray::Array<f64, ndarray::Ix2>>,
    pub alpha: RwLock<ndarray::Array<f64, ndarray::Ix2>>,
    pub beta: RwLock<ndarray::Array<f64, ndarray::Ix2>>,
    pub all_D: RwLock<ndarray::Array<f64, ndarray::Ix2>>,
    pub unique_guides: Vec<String>,
    pub unique_genes: Vec<String>,
    pub gene_to_guide_map: HashMap<String, Vec<String>>,
    pub guide_to_gene_map: HashMap<String, String>,
    pub lookup: RwLock<HashMap<&'a str, HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>>>>,
}

impl<'a> CAVI<'a> {
    pub fn new(
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
    ) -> Result<CAVI<'a>, Box<dyn std::error::Error>> {
        let all_x = RwLock::new(Array::from_shape_vec((all_x_nrow), all_x)?);
        let all_xx = RwLock::new(Array::from_shape_vec((all_x_nrow), all_xx)?);
        let all_y = RwLock::new(Array::from_shape_vec((all_y_nrow, all_y_ncol), all_y)?);
        let all_yy = RwLock::new(Array::from_shape_vec((all_y_nrow, all_y_ncol), all_yy)?);
        let all_r = RwLock::new(Array::from_shape_vec((all_r_nrow), all_r)?);
        let all_rr = RwLock::new(Array::from_shape_vec((all_r_nrow), all_rr)?);
        let all_s = RwLock::new(Array::from_shape_vec((all_s_nrow, all_s_ncol), all_s)?);
        let all_ss = RwLock::new(Array::from_shape_vec((all_s_nrow, all_s_ncol), all_ss)?);
        let all_tau = RwLock::new(Array::from_shape_vec(
            (all_tau_nrow, all_tau_ncol),
            all_tau,
        )?);
        let alpha = RwLock::new(Array::from_shape_vec(
            (all_tau_nrow, all_tau_ncol),
            all_alpha,
        )?);
        let beta = RwLock::new(Array::from_shape_vec(
            (all_tau_nrow, all_tau_ncol),
            all_beta,
        )?);
        let all_D = RwLock::new(Array::from_shape_vec((all_tau_nrow, all_tau_ncol), all_D)?);
        let MAE = RwLock::new(MAE);
        let lookup = RwLock::new(HashMap::with_capacity(6));
        let myCAVI = CAVI {
            verbose,
            debug,
            force,
            MAX_ITERATIONS,
            MAE_THRESHOLD,
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
            ncgenes,
            ncguides,
            sample_names,
            MAE,
            hash_x,
            all_x,
            all_xx,
            hash_y,
            all_y,
            all_yy,
            hash_r,
            all_r,
            all_rr,
            hash_s,
            all_s,
            all_ss,
            hash_td,
            all_tau,
            alpha,
            beta,
            all_D,
            unique_guides,
            unique_genes,
            gene_to_guide_map,
            guide_to_gene_map,
            lookup,
        };
        return Ok(myCAVI);
    }

    pub fn build_lookup(&'a self) -> Result<(), Box<dyn std::error::Error>> {
        if self.verbose {
            println!("Extracting references for speed up.");
            println!("Extracting x values.")
        }
        let x_value: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            self.generate_x_lookup_value();
        if self.verbose {println!("Extracting y values.");}
        let y_value: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            self.generate_y_lookup_value();
        if self.verbose {println!("Extracting r values.");}
        let r_value: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            self.generate_r_lookup_value();
        if self.verbose {println!("Extracting s values.");}
        let s_value: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            self.generate_s_lookup_value();
        if self.verbose {println!("Extracting tau values.");}
        let tau_value: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            self.generate_tau_lookup_value();
        if self.verbose {println!("Extracting mae values.");}
        let mae_value: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            self.generate_mae_lookup_value();
        self.lookup.write().unwrap().insert("x", x_value);
        self.lookup.write().unwrap().insert("y", y_value);
        self.lookup.write().unwrap().insert("r",r_value);
        self.lookup.write().unwrap().insert("s",s_value);
        self.lookup.write().unwrap().insert("tau",tau_value);
        self.lookup.write().unwrap().insert("mae",mae_value);
        return Ok(());
    }

    fn generate_x_lookup_value(&'a self) -> HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> {
        let mut out: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            HashMap::with_capacity(self.hash_x.keys().len());
        for gi in self.hash_x.keys() {
            let indices: Vec<HashMap<&'a str, &'a str>> = self.collect_all_x_keys(gi);
            out.insert(gi, indices);
        }
        return out;
    }

    fn generate_y_lookup_value(&'a self) -> HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> {
        let mut out: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            HashMap::with_capacity(self.hash_y.keys().len());
        for gi in self.hash_y.keys() {
            let indices: Vec<HashMap<&'a str, &'a str>> = self.collect_all_y_keys(gi);
            out.insert(gi, indices);
        }
        return out;
    }

    fn generate_r_lookup_value(&'a self) -> HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> {
        let mut out: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            HashMap::with_capacity(self.hash_r.keys().len());
        for gi in self.hash_r.keys() {
            let indices: Vec<HashMap<&'a str, &'a str>> = self.collect_all_r_keys(gi);
            out.insert(gi, indices);
        }
        return out;
    }

    fn generate_s_lookup_value(&'a self) -> HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> {
        let mut out: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            HashMap::with_capacity(self.hash_s.keys().len());
        for gi in self.hash_s.keys() {
            let indices: Vec<HashMap<&'a str, &'a str>> = self.collect_all_s_keys(gi);
            out.insert(gi, indices);
        }
        return out;
    }

    fn generate_tau_lookup_value(&'a self) -> HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> {
        let mut out: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            HashMap::with_capacity(self.hash_td.keys().len());
        for gi in self.hash_td.keys() {
            let indices: Vec<HashMap<&'a str, &'a str>> = self.collect_all_tau_keys(gi);
            out.insert(gi, indices);
        }
        return out;
    }

    fn generate_mae_lookup_value(&'a self) -> HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> {
        let mut out: HashMap<&'a str, Vec<HashMap<&'a str, &'a str>>> =
            HashMap::with_capacity(1);
        let indices: Vec<HashMap<&'a str, &'a str>> = self.collect_all_mae_keys();
        out.insert("all", indices);
        return out;
    }

    pub fn start_CAVI(&'a self) -> Result<(), Box<dyn std::error::Error>> {
        self.build_lookup()?;
        let mut counter = 0;
        while counter <= self.MAX_ITERATIONS {
            self.update_x()?;
            self.update_xx()?;
            self.update_y()?;
            self.update_yy()?;
            self.update_r()?;
            self.update_rr()?;
            self.update_s()?;
            self.update_ss()?;
            self.update_tau()?;
            self.update_MAE()?;
            counter = counter + 1;
            if self.verbose {
                println!("Iteration {} complete.", counter);
                println!(
                    "MAE at {}.",
                    self.MAE.read().unwrap()[self.MAE.read().unwrap().len() - 1]
                );
            }
            if self.check_convergence() {
                break;
            }
        }
        return Ok(());
    }

    fn check_convergence(&self) -> bool {
        let MAE_len = self.MAE.read().unwrap().len();
        if !self.force {
            if self.MAE.read().unwrap()[MAE_len - 1] > self.MAE.read().unwrap()[MAE_len - 2] {
                println!("Breaking CAVI: MAE increased with step. \n Init model with force=True if you want to keep going.");
                return true;
            }
        }
        if (self.MAE.read().unwrap()[MAE_len - 1] - self.MAE.read().unwrap()[MAE_len - 2]).abs()
            <= self.MAE_THRESHOLD
        {
            println!("Breaking CAVI: MAE_THRESHOLD reached.");
            return true;
        }
        if self.MAE.read().unwrap()[MAE_len - 1] < 0.0 {
            println!("Breaking CAVI: MAE turned negative.");
            return true;
        }
        return false;
    }
}
