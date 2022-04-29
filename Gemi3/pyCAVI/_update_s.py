import numpy as np
import sys
from Gemi3.pyCAVI.util import find_s_keys,find_r_keys,gene_to_guides_fast

def s_loop(self):
    if self.verbose:
        print("Updating s values.",file=sys.stdout)
    for s_triple in self.hash_s.keys():
        indices=collect_all_keys(self,s_triple)
        numerator_sum=np.sum(np.asarray(list(map(lambda x: s_numerator(self,x),indices))))
        denominator_sum=np.sum(np.asarray(list(map(lambda x: s_denominator(self,x),indices))))
        new_s=((self.prio_tau_s*self.prio_mu_s)+numerator_sum)/(self.prio_tau_s+denominator_sum)
        if self.debug: print(f"Old s {self.s[self.hash_s[s_triple]]}.", sys.stdout)
        self.s[self.hash_s[s_triple]] = new_s
        if self.debug: print(f"New s {self.s[self.hash_s[s_triple]]}.", sys.stdout)

def ss_loop(self):
    if self.verbose:
        print("Updating ss values.",file=sys.stdout)
    for s_triple in self.hash_s.keys():
        indices=collect_all_keys(self,s_triple)
        numerator_sum=np.sum(np.asarray(list(map(lambda x: s_numerator(self,x),indices))))
        denominator_sum=np.sum(np.asarray(list(map(lambda x: s_denominator(self,x),indices))))
        new_ss=((self.prio_tau_ss*self.prio_mu_ss)+numerator_sum)/(self.prio_tau_ss+denominator_sum)
        if self.debug: print(f"Old ss {self.ss[self.hash_s[s_triple]]}.", sys.stdout)
        self.ss[self.hash_s[s_triple]] = new_ss
        if self.debug: print(f"New ss {self.ss[self.hash_s[s_triple]]}.", sys.stdout)

def s_numerator(self,key_dict):
    numerator=(self.tau[self.hash_t[key_dict["t"]]]*self.r[self.hash_r[key_dict["r"]]])*\
    (self.D[self.hash_d[key_dict["t"]]]-\
    (self.x[self.hash_x[key_dict["hj"]]]*self.y[self.hash_y[key_dict["hl"]]])-\
    (self.x[self.hash_x[key_dict["fk"]]]*self.y[self.hash_y[key_dict["fl"]]])-\
    (self.x[self.hash_x[key_dict["gi"]]]*self.y[self.hash_y[key_dict["gl"]]]))
    return numerator

def s_denominator(self,key_dict):
    denominator = self.tau[self.hash_t[key_dict["t"]]] * self.rr[self.hash_r[key_dict["r"]]]
    return denominator

def collect_all_keys(self,s_triple):
    #only get keys from tau and split them accordingly
    gi, hj, fk, l = s_triple.split("_")
    s_triple="_".join([gi, hj, fk])
    s_as_guides=[gene_to_guides_fast(self,x) for x in s_triple.split("_")]
    s_as_guides = [x for y in s_as_guides for x in y]
    tau_keys = []
    for k in self.hash_t.keys():
        gi, hj, fk, l = k.split("_")
        if gi in s_as_guides and hj in s_as_guides and fk in s_as_guides:
            tau_keys.append(k)
        else:
            continue
    out=[]
    for k in tau_keys:
        gi,hj,fk,l = k.split("_")
        gl = self.guide_to_gene(gi)
        hl = self.guide_to_gene(hj)
        fl = self.guide_to_gene(fk)
        s = find_s_keys(self,gl,hl,fl,l)
        r = find_r_keys(self,gi,hj,fk)
        out.append({"t":k,"gi":gi,"hj":hj,"fk":fk,"l":l,
                    "gl":"_".join([gl,l]),
                    "hl":"_".join([hl,l]),
                    "fl":"_".join([hl,l]),
                    "r":r,
                    "s":s})
    return out