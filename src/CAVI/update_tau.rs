use super::CAVI;
use std::collections::HashMap;
use itertools::Itertools;
use rayon::prelude::*;

impl<'a> CAVI<'a>{
    pub fn update_tau(&self) -> Result<(),Box<dyn std::error::Error>>{
        if self.verbose {
            println! {"Updating tau values"}
        }
        for gl in self.lookup.read().unwrap().get("tau").unwrap().keys(){
            if self.lookup.read().unwrap().get("tau").unwrap().get(gl).unwrap().len() != 1 {panic!("Wrong length of indices in update_tau.rs")}
            let new_alpha = tau_alpha(self,&self.lookup.read().unwrap().get("tau").unwrap().get(gl).unwrap()[0]);
            {*self.alpha.write().unwrap().get_mut(*self.hash_td.get(*gl).unwrap()).unwrap()= new_alpha}
            let new_beta = tau_beta(self,&self.lookup.read().unwrap().get("tau").unwrap().get(gl).unwrap()[0]);
            {*self.beta.write().unwrap().get_mut(*self.hash_td.get(*gl).unwrap()).unwrap()= new_beta}
            let new_tau = *self.alpha.read().unwrap().get(*self.hash_td.get(*gl).unwrap()).unwrap() / *self.beta.read().unwrap().get(*self.hash_td.get(*gl).unwrap()).unwrap();
            if self.debug {println!("Old tau {}.",*self.all_tau.read().unwrap().get(*self.hash_td.get(*gl).unwrap()).unwrap())}
            {*self.all_tau.write().unwrap().get_mut(*self.hash_td.get(*gl).unwrap()).unwrap()= new_tau}
            if self.debug {println!("New tau {}.",*self.all_tau.read().unwrap().get(*self.hash_td.get(*gl).unwrap()).unwrap())}
        }
        return Ok(());
    }
}


fn tau_alpha(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let out : f64 = *obj.alpha.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() + 0.5;
    return out;
    }

fn tau_beta(obj : &CAVI,key_dict : &HashMap<&str,&str>) -> f64 {
    let beta : f64 = *obj.beta.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap();
    let D2 : f64 = (0.5 * *obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap()).powf(2.0);
    let Dxgiygl : f64 = *obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() * *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap()).unwrap();
    let Dxhjyhl : f64 = *obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() * *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap();
    let Dxfkyfl : f64 = *obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() * *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("fk").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("fl").unwrap()).unwrap()).unwrap();
    let Drs : f64 = *obj.all_D.read().unwrap().get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap()).unwrap() * *obj.all_r.read().unwrap().get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap()).unwrap() * *obj.all_s.read().unwrap().get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap()).unwrap();
    let xgi2ygl2 :f64 = *obj.all_xx.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap() * *obj.all_yy.read().unwrap().get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap()).unwrap();
    let xgiyglxhjyhl :f64 = *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap()).unwrap() * *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap();
    let xgiyglxfkyfl : f64 =*obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap()).unwrap() * *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("fk").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("fl").unwrap()).unwrap()).unwrap();
    let xgiyglrs :f64 = 2.0 * *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("gi").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap()).unwrap() * *obj.all_r.read().unwrap().get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap()).unwrap() * *obj.all_s.read().unwrap().get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap()).unwrap();
    let xhj2yhl2 :f64 = 0.5 * *obj.all_xx.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_yy.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap();
    let xhjyhlxfkyfl :f64 = *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap() * *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("fk").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("fl").unwrap()).unwrap()).unwrap();
    let xhjyhlrs :f64 = *obj.all_x.read().unwrap().get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap()).unwrap() * *obj.all_y.read().unwrap().get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap()).unwrap() * *obj.all_r.read().unwrap().get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap()).unwrap() * *obj.all_s.read().unwrap().get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap()).unwrap();
    let r2s2 :f64 = 0.5 * *obj.all_rr.read().unwrap().get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap()).unwrap() * *obj.all_ss.read().unwrap().get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap()).unwrap();
    let out : f64 = beta+D2-Dxgiygl-Dxhjyhl-Dxfkyfl-Drs+xgi2ygl2+xgiyglxhjyhl+xgiyglxfkyfl+xgiyglrs+xhj2yhl2+xhjyhlxfkyfl+xhjyhlrs+r2s2;
    return out;
    }
