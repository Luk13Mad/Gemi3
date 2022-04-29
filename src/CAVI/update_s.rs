use super::CAVI;
use std::collections::HashMap;
use itertools::Itertools;
use rayon::prelude::*;

impl<'a> CAVI<'a>{
    pub fn update_s(&self) -> Result<(),Box<dyn std::error::Error>>{
        if self.verbose {
            println! {"Updating s values"}
        }
        for gl in self.lookup.read().unwrap().get("s").unwrap().keys(){
            let all_numerators :Vec<f64> = self.lookup.read().unwrap().get("s").unwrap().get(gl).unwrap().par_iter().map(|element| s_numerator(self,element)).collect();
            let all_denominators :Vec<f64> = self.lookup.read().unwrap().get("s").unwrap().get(gl).unwrap().par_iter().map(|element| s_denominator(self,element)).collect();
            let new_s : f64 = ((self.prio_tau_s * self.prio_mu_s) + all_numerators.par_iter().sum::<f64>())/(self.prio_tau_s+all_denominators.par_iter().sum::<f64>());
            if self.debug {println!("Old s {}.",*self.all_s.read().unwrap().get(*self.hash_s.get(*gl).unwrap()).unwrap())}
            {*self.all_s.write().unwrap().get_mut(*self.hash_s.get(*gl).unwrap()).unwrap()= new_s}
            if self.debug {println!("New s {}.",*self.all_s.read().unwrap().get(*self.hash_s.get(*gl).unwrap()).unwrap())}
        }
        return Ok(());
    }

    pub fn update_ss(&self) -> Result<(),Box<dyn std::error::Error>>{
        if self.verbose{
            println!{"Updating ss values"}
        }
        for gl in self.lookup.read().unwrap().get("s").unwrap().keys(){
            let all_numerators :Vec<f64> = self.lookup.read().unwrap().get("s").unwrap().get(gl).unwrap().par_iter().map(|element| s_numerator(self,element)).collect();
            let all_denominators :Vec<f64> = self.lookup.read().unwrap().get("s").unwrap().get(gl).unwrap().par_iter().map(|element| s_denominator(self,element)).collect();
            let new_ss : f64 = ((self.prio_tau_ss * self.prio_mu_ss) + all_numerators.par_iter().sum::<f64>())/(self.prio_tau_ss+all_denominators.par_iter().sum::<f64>());
            if self.debug {println!("Old ss {}.",*self.all_ss.read().unwrap().get(*self.hash_s.get(*gl).unwrap()).unwrap())}
            {*self.all_ss.write().unwrap().get_mut(*self.hash_s.get(*gl).unwrap()).unwrap()= new_ss}
            if self.debug {println!("New ss {}.",*self.all_ss.read().unwrap().get(*self.hash_s.get(*gl).unwrap()).unwrap())}
        }
        return Ok(());
    }

}

fn s_numerator(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let numerator : f64 = (*obj.all_tau.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() *
        *obj.all_r.read().unwrap().get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap()).unwrap()) *
        (*obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("fk").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("fl").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap()).unwrap());
    return numerator;
    }

fn s_denominator(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let denominator : f64 = *obj.all_tau.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() *
        *obj.all_rr.read().unwrap().get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap()).unwrap();
    return denominator;
    }
