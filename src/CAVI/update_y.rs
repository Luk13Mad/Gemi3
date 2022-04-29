use super::CAVI;
use std::collections::HashMap;
use itertools::Itertools;
use rayon::prelude::*;

impl<'a> CAVI<'a>{
    pub fn update_y(&self) -> Result<(),Box<dyn std::error::Error>>{
        if self.verbose {
            println! {"Updating y values"}
        }
        for gl in self.lookup.read().unwrap().get("y").unwrap().keys(){
            let all_numerators :Vec<f64> = self.lookup.read().unwrap().get("y").unwrap().get(gl).unwrap().par_iter().map(|element| y_numerator(self,element)).collect();
            let all_denominators :Vec<f64> = self.lookup.read().unwrap().get("y").unwrap().get(gl).unwrap().par_iter().map(|element| y_denominator(self,element)).collect();
            let new_ygl : f64 = ((self.prio_tau_y * self.prio_mu_y) + all_numerators.par_iter().sum::<f64>())/(self.prio_tau_y+all_denominators.par_iter().sum::<f64>());
            if self.debug {println!("Old ygl {}.",*self.all_y.read().unwrap().get(*self.hash_y.get(*gl).unwrap()).unwrap())}
            {*self.all_y.write().unwrap().get_mut(*self.hash_y.get(*gl).unwrap()).unwrap()= new_ygl}
            if self.debug {println!("New ygl {}.",*self.all_y.read().unwrap().get(*self.hash_y.get(*gl).unwrap()).unwrap())}
        }
        return Ok(());
    }

    pub fn update_yy(&self) -> Result<(),Box<dyn std::error::Error>>{
        if self.verbose{
            println!{"Updating yy values"}
        }
        for gl in self.lookup.read().unwrap().get("y").unwrap().keys(){
            let all_numerators :Vec<f64> = self.lookup.read().unwrap().get("y").unwrap().get(gl).unwrap().par_iter().map(|element| y_numerator(self,element)).collect();
            let all_denominators :Vec<f64> = self.lookup.read().unwrap().get("y").unwrap().get(gl).unwrap().par_iter().map(|element| y_denominator(self,element)).collect();
            let new_yygl : f64 = ((self.prio_tau_yy * self.prio_mu_yy) + all_numerators.par_iter().sum::<f64>())/(self.prio_tau_yy+all_denominators.par_iter().sum::<f64>());
            if self.debug {println!("Old yygl {}.",*self.all_yy.read().unwrap().get(*self.hash_y.get(*gl).unwrap()).unwrap())}
            {*self.all_yy.write().unwrap().get_mut(*self.hash_y.get(*gl).unwrap()).unwrap()= new_yygl}
            if self.debug {println!("New yygl {}.",*self.all_yy.read().unwrap().get(*self.hash_y.get(*gl).unwrap()).unwrap())}
        }
        return Ok(());
    }

}

fn y_numerator(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let numerator : f64 = (*obj.all_tau.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() *
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap()) *
        (*obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap() -
        *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("fk").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("fl").unwrap()).unwrap()).unwrap() -
        *obj.all_r.read().unwrap().get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap()).unwrap() * *obj.all_s.read().unwrap().get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap()).unwrap());
    return numerator;
    }


fn y_denominator(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let denominator : f64 = *obj.all_tau.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() *
        *obj.all_xx.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap();
    return denominator;
    }
