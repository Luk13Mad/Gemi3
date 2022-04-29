use super::CAVI;
use rayon::prelude::*;
use std::collections::HashMap;

impl<'a> CAVI<'a> {
    pub fn update_x(&self) -> Result<(), Box<dyn std::error::Error>> {
        if self.verbose {
            println! {"Updating x values"}
        }
        for gi in self.lookup.read().unwrap().get("x").unwrap().keys() {
            let all_numerators: Vec<f64> = self
                .lookup
                .read()
                .unwrap()
                .get("x")
                .unwrap()
                .get(gi)
                .unwrap()
                .par_iter()
                .map(|element| x_numerator(self, element))
                .collect();
            let all_denominators: Vec<f64> = self
                .lookup
                .read()
                .unwrap()
                .get("x")
                .unwrap()
                .get(gi)
                .unwrap()
                .par_iter()
                .map(|element| x_denominator(self, element))
                .collect();
            let new_xgi: f64 = ((self.prio_tau_x * self.prio_mu_x)
                + all_numerators.par_iter().sum::<f64>())
                / (self.prio_tau_x + all_denominators.par_iter().sum::<f64>());
            if self.debug {
                println!(
                    "Old xgi {}.",
                    *self
                        .all_x
                        .read()
                        .unwrap()
                        .get(*self.hash_x.get(*gi).unwrap())
                        .unwrap()
                )
            }
            {
                *self
                    .all_x
                    .write()
                    .unwrap()
                    .get_mut(*self.hash_x.get(*gi).unwrap())
                    .unwrap() = new_xgi
            }
            if self.debug {
                println!(
                    "New xgi {}.",
                    *self
                        .all_x
                        .read()
                        .unwrap()
                        .get(*self.hash_x.get(*gi).unwrap())
                        .unwrap()
                )
            }
        }
        return Ok(());
    }

    pub fn update_xx(&self) -> Result<(), Box<dyn std::error::Error>> {
        if self.verbose {
            println! {"Updating xx values"}
        }
        for gi in self.lookup.read().unwrap().get("x").unwrap().keys() {
            let all_numerators: Vec<f64> = self
                .lookup
                .read()
                .unwrap()
                .get("x")
                .unwrap()
                .get(gi)
                .unwrap()
                .par_iter()
                .map(|element| x_numerator(self, element))
                .collect();
            let all_denominators: Vec<f64> = self
                .lookup
                .read()
                .unwrap()
                .get("x")
                .unwrap()
                .get(gi)
                .unwrap()
                .par_iter()
                .map(|element| x_denominator(self, element))
                .collect();
            let new_xxgi: f64 = ((self.prio_tau_xx * self.prio_mu_xx)
                + all_numerators.par_iter().sum::<f64>())
                / (self.prio_tau_xx + all_denominators.par_iter().sum::<f64>());
            if self.debug {
                println!(
                    "Old xxgi {}.",
                    *self
                        .all_xx
                        .read()
                        .unwrap()
                        .get(*self.hash_x.get(*gi).unwrap())
                        .unwrap()
                )
            }
            {
                *self
                    .all_xx
                    .write()
                    .unwrap()
                    .get_mut(*self.hash_x.get(*gi).unwrap())
                    .unwrap() = new_xxgi
            }
            if self.debug {
                println!(
                    "New xxgi {}.",
                    *self
                        .all_xx
                        .read()
                        .unwrap()
                        .get(*self.hash_x.get(*gi).unwrap())
                        .unwrap()
                )
            }
        }
        return Ok(());
    }
}

fn x_numerator(obj: &CAVI, key_dict: &HashMap<&str, &str>) -> f64 {
    let numerator: f64 = (*obj
        .all_tau
        .read()
        .unwrap()
        .get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap())
        .unwrap()
        * *obj
            .all_y
            .read()
            .unwrap()
            .get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap())
            .unwrap())
        * (*obj
            .all_D
            .read()
            .unwrap()
            .get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap())
            .unwrap()
            - *obj
                .all_x
                .read()
                .unwrap()
                .get(*obj.hash_x.get(*key_dict.get("hj").unwrap()).unwrap())
                .unwrap()
                * *obj
                    .all_y
                    .read()
                    .unwrap()
                    .get(*obj.hash_y.get(*key_dict.get("hl").unwrap()).unwrap())
                    .unwrap()
            - *obj
                .all_x
                .read()
                .unwrap()
                .get(*obj.hash_x.get(*key_dict.get("fk").unwrap()).unwrap())
                .unwrap()
                * *obj
                    .all_y
                    .read()
                    .unwrap()
                    .get(*obj.hash_y.get(*key_dict.get("fl").unwrap()).unwrap())
                    .unwrap()
            - *obj
                .all_r
                .read()
                .unwrap()
                .get(*obj.hash_r.get(*key_dict.get("r").unwrap()).unwrap())
                .unwrap()
                * *obj
                    .all_s
                    .read()
                    .unwrap()
                    .get(*obj.hash_s.get(*key_dict.get("s").unwrap()).unwrap())
                    .unwrap());
    return numerator;
}

fn x_denominator(obj: &CAVI, key_dict: &HashMap<&str, &str>) -> f64 {
    let denominator: f64 = *obj
        .all_tau
        .read()
        .unwrap()
        .get(*obj.hash_td.get(*key_dict.get("t").unwrap()).unwrap())
        .unwrap()
        * *obj
            .all_yy
            .read()
            .unwrap()
            .get(*obj.hash_y.get(*key_dict.get("gl").unwrap()).unwrap())
            .unwrap();
    return denominator;
}
