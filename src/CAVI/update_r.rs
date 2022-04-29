use super::CAVI;
use std::collections::HashMap;
use itertools::Itertools;
use rayon::prelude::*;

impl<'a> CAVI<'a>{
    pub fn update_r(&self) -> Result<(),Box<dyn std::error::Error>>{
        if self.verbose {
            println! {"Updating r values"}
        }
        for gl in self.lookup.read().unwrap().get("r").unwrap().keys(){
            let all_numerators :Vec<f64> = self.lookup.read().unwrap().get("r").unwrap().get(gl).unwrap().par_iter().map(|element| r_numerator(self,element)).collect();
            let all_denominators :Vec<f64> = self.lookup.read().unwrap().get("r").unwrap().get(gl).unwrap().par_iter().map(|element| r_denominator(self,element)).collect();
            let new_r : f64 = ((self.prio_tau_r * self.prio_mu_r) + all_numerators.par_iter().sum::<f64>())/(self.prio_tau_r+all_denominators.par_iter().sum::<f64>());
            if self.debug {println!("Old r {}.",*self.all_r.read().unwrap().get(*self.hash_r.get(*gl).unwrap()).unwrap())}
            {*self.all_r.write().unwrap().get_mut(*self.hash_r.get(*gl).unwrap()).unwrap()= new_r}
            if self.debug {println!("New r {}.",*self.all_r.read().unwrap().get(*self.hash_r.get(*gl).unwrap()).unwrap())}
        }
        return Ok(());
    }

    pub fn update_rr(&self) -> Result<(),Box<dyn std::error::Error>>{
        if self.verbose{
            println!{"Updating rr values"}
        }
        for gl in self.lookup.read().unwrap().get("r").unwrap().keys(){
            let all_numerators :Vec<f64> = self.lookup.read().unwrap().get("r").unwrap().get(gl).unwrap().par_iter().map(|element| r_numerator(self,element)).collect();
            let all_denominators :Vec<f64> = self.lookup.read().unwrap().get("r").unwrap().get(gl).unwrap().par_iter().map(|element| r_denominator(self,element)).collect();
            let new_rr : f64 = ((self.prio_tau_rr * self.prio_mu_rr) + all_numerators.par_iter().sum::<f64>())/(self.prio_tau_rr+all_denominators.par_iter().sum::<f64>());
            if self.debug {println!("Old rr {}.",*self.all_rr.read().unwrap().get(*self.hash_r.get(*gl).unwrap()).unwrap())}
            {*self.all_rr.write().unwrap().get_mut(*self.hash_r.get(*gl).unwrap()).unwrap()= new_rr}
            if self.debug {println!("New rr {}.",*self.all_rr.read().unwrap().get(*self.hash_r.get(*gl).unwrap()).unwrap())}
        }
        return Ok(());
    }

}

fn r_numerator(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let numerator : f64 = (*obj.all_tau.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() *
        *obj.all_s.read().unwrap().get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap()).unwrap()) *
        (*obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("fk").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("fl").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap()).unwrap());
    return numerator;
    }


fn r_denominator(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let denominator : f64 = *obj.all_tau.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() *
        *obj.all_ss.read().unwrap().get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap()).unwrap();
    return denominator;
    }