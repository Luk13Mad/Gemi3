use super::CAVI;
use itertools::Itertools;

impl<'a> CAVI<'a> {
    pub fn guide_to_gene(&self, guide: &str) -> &str {
        let out = match self.guide_to_gene_map.get(guide) {
            Some(v) => v,
            None => panic!("Failed to map guide to gene"),
        };
        return &*out;
    }

    pub fn gene_to_guides(&self, gene: &str) -> &Vec<String> {
        let out = match self.gene_to_guide_map.get(gene) {
            Some(v) => v,
            None => panic!("Failed to map gene to guides"),
        };
        return &*out;
    }

    pub fn find_s_keys(&self, gl: &str, hl: &str, fl: &str, l: &str) -> Option<&str> {
        let concat = [gl, hl, fl, l].join("_");
        if self.hash_s.contains_key(&concat) {
            let pair = self.hash_s.get_key_value(&concat).unwrap();
            return Some(&*pair.0);
        } else {
            for k in self.hash_s.keys() {
                let (first, second, third, sample) = match k.split("_").collect_tuple() {
                    Some(v) => v,
                    None => panic!("Error find_s_key in utility.rs"),
                };
                if sample != l {
                    continue;
                } else if gl == first && hl == second && fl == third {
                    //123
                    return Some(k);
                } else if gl == second && hl == first && fl == third {
                    //213
                    return Some(k);
                } else if gl == first && hl == third && fl == second {
                    //132
                    return Some(k);
                } else if gl == second && hl == third && fl == first {
                    //231
                    return Some(k);
                } else if gl == third && hl == first && fl == second {
                    //312
                    return Some(k);
                } else if gl == third && hl == second && fl == first {
                    //321
                    return Some(k);
                }
            }
            return None;
        }
    }

    pub fn find_r_keys(&self, gi: &str, hj: &str, fk: &str) -> Option<&str> {
        let concat = [gi, hj, fk].join("_");
        if self.hash_r.contains_key(&concat) {
            let pair = self.hash_r.get_key_value(&concat).unwrap();
            return Some(&*pair.0);
        } else {
            for k in self.hash_r.keys() {
                let (first, second, third) = match k.split("_").collect_tuple() {
                    Some(v) => v,
                    None => panic!("Error find_r_key in utility.rs"),
                };
                if gi == first && hj == second && fk == third {
                    //123
                    return Some(k);
                } else if gi == second && hj == first && fk == third {
                    //213
                    return Some(k);
                } else if gi == first && hj == third && fk == second {
                    //132
                    return Some(k);
                } else if gi == second && hj == third && fk == first {
                    //231
                    return Some(k);
                } else if gi == third && hj == first && fk == second {
                    //312
                    return Some(k);
                } else if gi == third && hj == second && fk == first {
                    //321
                    return Some(k);
                }
            }
            return None;
        }
    }
}
