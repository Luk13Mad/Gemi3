use super::CAVI;
use itertools::Itertools;
use std::collections::HashMap;

impl<'a> CAVI<'a> {
    pub fn collect_all_x_keys(&'a self, key: &'a str) -> Vec<HashMap<&'a str, &'a str>> {
        let mut out: Vec<HashMap<&str, &str>> = Vec::new();
        let tau_keys: Vec<&'a str> = self
            .hash_td
            .keys()
            .filter(|element| element.contains(key))
            .map(|element| &element[..])
            .collect();
        if self.debug {
            println!("update x : {:?}", &tau_keys)
        }
        for k in tau_keys {
            let (gi, hj, fk, l) = match k.split("_").collect_tuple() {
                Some(v) => v,
                None => panic!("Error collect_all_x_keys"),
            };
            let gl = self.guide_to_gene(&gi);
            let hl = self.guide_to_gene(&hj);
            let fl = self.guide_to_gene(&fk);
            let s = match self.find_s_keys(&gl, &hl, &fl, &l) {
                Some(v) => v,
                None => panic!("Error find_s_keys in collect_all_x_keys"),
            };
            let r = match self.find_r_keys(&gi, &hj, &fk) {
                Some(v) => v,
                None => panic!("Error find_r_keys in collect_all_x_keys"),
            };
            let mut internal_hashmap: HashMap<&str, &str> = HashMap::with_capacity(10);
            internal_hashmap.insert("t", &k);
            internal_hashmap.insert("gi", &gi);
            internal_hashmap.insert("hj", &hj);
            internal_hashmap.insert("fk", &fk);
            internal_hashmap.insert("l", l);
            internal_hashmap.insert(
                "gl",
                &*(self.hash_y.get_key_value(&[gl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "hl",
                &*(self.hash_y.get_key_value(&[hl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "fl",
                &*(self.hash_y.get_key_value(&[fl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert("r", r);
            internal_hashmap.insert("s", s);
            out.push(internal_hashmap);
        }
        return out;
    }

    pub fn collect_all_y_keys(&'a self, key: &'a str) -> Vec<HashMap<&'a str, &'a str>> {
        let mut out: Vec<HashMap<&str, &str>> = Vec::new();
        let (gene,l) = match key.split("_").collect_tuple(){Some(v) => v, None => panic!("Error collect_all_y_keys.")};
        let key_guides : &Vec<String> = self.gene_to_guides(&gene);
        let mut tau_keys : Vec<&'a str> = Vec::new();
        for k in key_guides{
            tau_keys.append(&mut self.hash_td.keys().filter(|element| element.contains(k)).map(|element| &element[..]).collect());
        }
        if self.debug {
            println!("update y : {:?}", &tau_keys)
        }
        for k in tau_keys {
            let (gi, hj, fk, l) = match k.split("_").collect_tuple() {
                Some(v) => v,
                None => panic!("Error collect_all_keys in collect_all_y_keys"),
            };
            let gl = self.guide_to_gene(&gi);
            let hl = self.guide_to_gene(&hj);
            let fl = self.guide_to_gene(&fk);
            let s = match self.find_s_keys(&gl, &hl, &fl, &l) {
                Some(v) => v,
                None => panic!("Error find_s_keys in collect_all_y_keys"),
            };
            let r = match self.find_r_keys(&gi, &hj, &fk) {
                Some(v) => v,
                None => panic!("Error find_r_keys in collect_all_y_keys"),
            };
            let mut internal_hashmap: HashMap<&str, &str> = HashMap::with_capacity(10);
            internal_hashmap.insert("t", &k);
            internal_hashmap.insert("gi", &gi);
            internal_hashmap.insert("hj", &hj);
            internal_hashmap.insert("fk", &fk);
            internal_hashmap.insert("l", l);
            internal_hashmap.insert(
                "gl",
                &*(self.hash_y.get_key_value(&[gl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "hl",
                &*(self.hash_y.get_key_value(&[hl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "fl",
                &*(self.hash_y.get_key_value(&[fl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert("r", r);
            internal_hashmap.insert("s", s);
            out.push(internal_hashmap);
        }
        return out;
    }

    pub fn collect_all_r_keys(&'a self, key: &'a str) -> Vec<HashMap<&'a str, &'a str>> {
        let mut out: Vec<HashMap<&str, &str>> = Vec::new();
        let (first_r,second_r,third_r) = match key.split("_").collect_tuple(){Some(v) => v, None => panic!("Error collect_all_r_keys")};
        //let tau_keys : Vec<String> = obj.hash_td.keys().filter(|element| element.contains(key)).map(|element| element.clone()).collect();
        let mut key_guides : Vec<&str> = Vec::with_capacity(3);
        key_guides.push(first_r);
        key_guides.push(second_r);
        key_guides.push(third_r);

        let mut tau_keys : Vec<&'a str> = Vec::new();
        for k in self.hash_td.keys() {
            let (first_t, second_t, third_t, l) = match k.split("_").collect_tuple() {
                Some(v) => v,
                None => panic!("Error collect_all_r_keys")
            };
            if key_guides.contains(&first_t) && key_guides.contains(&second_t) && key_guides.contains(&third_t) {
                tau_keys.push(&*k);
            }
        }

        if self.debug {
            println!("update r : {:?}", &tau_keys)
        }
        for k in tau_keys {
            let (gi, hj, fk, l) = match k.split("_").collect_tuple() {
                Some(v) => v,
                None => panic!("Error collect_all_keys in collect_all_r_keys"),
            };
            let gl = self.guide_to_gene(&gi);
            let hl = self.guide_to_gene(&hj);
            let fl = self.guide_to_gene(&fk);
            let s = match self.find_s_keys(&gl, &hl, &fl, &l) {
                Some(v) => v,
                None => panic!("Error find_s_keys in collect_all_r_keys"),
            };
            let r = match self.find_r_keys(&gi, &hj, &fk) {
                Some(v) => v,
                None => panic!("Error find_r_keys in collect_all_r_keys"),
            };
            let mut internal_hashmap: HashMap<&str, &str> = HashMap::with_capacity(10);
            internal_hashmap.insert("t", &k);
            internal_hashmap.insert("gi", &gi);
            internal_hashmap.insert("hj", &hj);
            internal_hashmap.insert("fk", &fk);
            internal_hashmap.insert("l", l);
            internal_hashmap.insert(
                "gl",
                &*(self.hash_y.get_key_value(&[gl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "hl",
                &*(self.hash_y.get_key_value(&[hl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "fl",
                &*(self.hash_y.get_key_value(&[fl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert("r", r);
            internal_hashmap.insert("s", s);
            out.push(internal_hashmap);
        }
        return out;
    }

    pub fn collect_all_s_keys(&'a self, key: &'a str) -> Vec<HashMap<&'a str, &'a str>> {
        let mut out: Vec<HashMap<&str, &str>> = Vec::new();
        let (first_s,second_s,third_s,l_s) = match key.split("_").collect_tuple(){Some(v) => v, None => panic!("Error collect_all_s_keys")};
        let mut key_guides : Vec<String> = self.gene_to_guides(first_s).clone();
        key_guides.append(&mut self.gene_to_guides(second_s).clone());
        key_guides.append(&mut self.gene_to_guides(third_s).clone());

        let mut tau_keys : Vec<&'a str> = Vec::new();
        for k in self.hash_td.keys(){
            let (first_t,second_t,third_t,l_t) = match k.split("_").collect_tuple(){Some(v) => v, None => panic!("Error collect_all_s_keys")};
            if key_guides.contains(&first_t.to_string()) &&  key_guides.contains(&second_t.to_string()) && key_guides.contains(&third_t.to_string()) && (l_s == l_t) {
            tau_keys.push(&*k);
            }
        }

        if self.debug {
            println!("update s : {:?}", &tau_keys)
        }
        for k in tau_keys {
            let (gi, hj, fk, l) = match k.split("_").collect_tuple() {
                Some(v) => v,
                None => panic!("Error collect_all_s_keys"),
            };
            let gl = self.guide_to_gene(&gi);
            let hl = self.guide_to_gene(&hj);
            let fl = self.guide_to_gene(&fk);
            let s = match self.find_s_keys(&gl, &hl, &fl, &l) {
                Some(v) => v,
                None => panic!("Error find_s_keys in collect_all_s_keys"),
            };
            let r = match self.find_r_keys(&gi, &hj, &fk) {
                Some(v) => v,
                None => panic!("Error find_r_keys in collect_all_s_keys"),
            };
            let mut internal_hashmap: HashMap<&str, &str> = HashMap::with_capacity(10);
            internal_hashmap.insert("t", &k);
            internal_hashmap.insert("gi", &gi);
            internal_hashmap.insert("hj", &hj);
            internal_hashmap.insert("fk", &fk);
            internal_hashmap.insert("l", l);
            internal_hashmap.insert(
                "gl",
                &*(self.hash_y.get_key_value(&[gl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "hl",
                &*(self.hash_y.get_key_value(&[hl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "fl",
                &*(self.hash_y.get_key_value(&[fl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert("r", r);
            internal_hashmap.insert("s", s);
            out.push(internal_hashmap);
        }
        return out;
    }

    pub fn collect_all_tau_keys(&'a self, key: &'a str) -> Vec<HashMap<&'a str, &'a str>> {
        let mut out: Vec<HashMap<&str, &str>> = Vec::new();
        let tau_keys : Vec<&'a str> = vec![key];

        if self.debug {
            println!("update s : {:?}", &tau_keys)
        }
        for k in tau_keys {
            let (gi, hj, fk, l) = match k.split("_").collect_tuple() {
                Some(v) => v,
                None => panic!("Error collect_all_s_keys"),
            };
            let gl = self.guide_to_gene(&gi);
            let hl = self.guide_to_gene(&hj);
            let fl = self.guide_to_gene(&fk);
            let s = match self.find_s_keys(&gl, &hl, &fl, &l) {
                Some(v) => v,
                None => panic!("Error find_s_keys in collect_all_s_keys"),
            };
            let r = match self.find_r_keys(&gi, &hj, &fk) {
                Some(v) => v,
                None => panic!("Error find_r_keys in collect_all_s_keys"),
            };
            let mut internal_hashmap: HashMap<&str, &str> = HashMap::with_capacity(10);
            internal_hashmap.insert("t", &k);
            internal_hashmap.insert("gi", &gi);
            internal_hashmap.insert("hj", &hj);
            internal_hashmap.insert("fk", &fk);
            internal_hashmap.insert("l", l);
            internal_hashmap.insert(
                "gl",
                &*(self.hash_y.get_key_value(&[gl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "hl",
                &*(self.hash_y.get_key_value(&[hl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "fl",
                &*(self.hash_y.get_key_value(&[fl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert("r", r);
            internal_hashmap.insert("s", s);
            out.push(internal_hashmap);
        }
        return out;
    }

    pub fn collect_all_mae_keys(&'a self) -> Vec<HashMap<&'a str, &'a str>> {
        let mut out: Vec<HashMap<&str, &str>> = Vec::new();
        let tau_keys : Vec<&str> = self.hash_td.keys().map(|key| &key[..]).collect();
        if self.debug {
            println!("update s : {:?}", &tau_keys)
        }
        for k in tau_keys {
            let (gi, hj, fk, l) = match k.split("_").collect_tuple() {
                Some(v) => v,
                None => panic!("Error collect_all_s_keys"),
            };
            let gl = self.guide_to_gene(&gi);
            let hl = self.guide_to_gene(&hj);
            let fl = self.guide_to_gene(&fk);
            let s = match self.find_s_keys(&gl, &hl, &fl, &l) {
                Some(v) => v,
                None => panic!("Error find_s_keys in collect_all_s_keys"),
            };
            let r = match self.find_r_keys(&gi, &hj, &fk) {
                Some(v) => v,
                None => panic!("Error find_r_keys in collect_all_s_keys"),
            };
            let mut internal_hashmap: HashMap<&str, &str> = HashMap::with_capacity(10);
            internal_hashmap.insert("t", &k);
            internal_hashmap.insert("gi", &gi);
            internal_hashmap.insert("hj", &hj);
            internal_hashmap.insert("fk", &fk);
            internal_hashmap.insert("l", l);
            internal_hashmap.insert(
                "gl",
                &*(self.hash_y.get_key_value(&[gl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "hl",
                &*(self.hash_y.get_key_value(&[hl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert(
                "fl",
                &*(self.hash_y.get_key_value(&[fl, l].join("_")).unwrap().0),
            );
            internal_hashmap.insert("r", r);
            internal_hashmap.insert("s", s);
            out.push(internal_hashmap);
        }
        return out;
    }

}

