use super::CAVI;
use std::collections::HashMap;
use itertools::Itertools;
use rayon::prelude::*;

impl<'a> CAVI<'a>{
    pub fn update_MAE(&self) -> Result<(),Box<dyn std::error::Error>>{
        if self.verbose {
            println! {"Updating MAE value"}
        }
        let new_MAE_vec : Vec<f64> = self.lookup.read().unwrap().get("mae").unwrap().get("all").unwrap().par_iter().map(|element| MAE_one_triple(self,element)).collect();
        let new_MAE : f64 = calc_MAE(&new_MAE_vec);
        {self.MAE.write().unwrap().push(new_MAE)}
        if self.debug {println!("New MAE {}.",self.MAE.read().unwrap()[self.MAE.read().unwrap().len()-1])}
        return Ok(());
    }
}

fn calc_MAE(values : &Vec<f64>) -> f64 {
    let sum : f64 = values.par_iter().map(|element| (*element).abs() ).sum::<f64>();
    let mean : f64 = sum / values.len() as f64;
    return mean;
}

fn MAE_one_triple(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let out : f64 = *obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("fk").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("fl").unwrap()).unwrap()).unwrap() -
        *obj.all_r.read().unwrap().get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap()).unwrap() * *obj.all_s.read().unwrap().get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap()).unwrap();
    return out;
    }
