import numpy as np
import sys
from Gemi3.pyCAVI.util import find_s_keys,find_r_keys

def r_loop(self):
    if self.verbose:
        print("Updating r values",file=sys.stdout)
    for r_triple in self.hash_r.keys():
        indices=collect_all_keys(self,r_triple)
        numerator_sum=np.sum(np.asarray(list(map(lambda x: r_numerator(self,x),indices))))
        denominator_sum=np.sum(np.asarray(list(map(lambda x: r_denominator(self,x),indices))))
        new_r=((self.prio_tau_r*self.prio_mu_r)+numerator_sum)/(self.prio_tau_r+denominator_sum)
        if self.debug: print(f"Old r {self.r[self.hash_r[r_triple]]}.", sys.stdout)
        self.r[self.hash_r[r_triple]]= new_r
        if self.debug: print(f"New r {self.r[self.hash_r[r_triple]]}.", sys.stdout)

def rr_loop(self):
    if self.verbose:
        print("Updating rr values",file=sys.stdout)
    for r_triple in self.hash_r.keys():
        indices=collect_all_keys(self,r_triple)
        numerator_sum=np.sum(np.asarray(list(map(lambda x: r_numerator(self,x),indices))))
        denominator_sum=np.sum(np.asarray(list(map(lambda x: r_denominator(self,x),indices))))
        new_rr=((self.prio_tau_rr*self.prio_mu_rr)+numerator_sum)/(self.prio_tau_rr+denominator_sum)
        if self.debug: print(f"Old rr {self.rr[self.hash_r[r_triple]]}.", sys.stdout)
        self.rr[self.hash_r[r_triple]]= new_rr
        if self.debug: print(f"New rr {self.rr[self.hash_r[r_triple]]}.", sys.stdout)

def r_numerator(self,key_dict):
    numerator = (self.tau[self.hash_t[key_dict["t"]]] * self.s[self.hash_s[key_dict["s"]]]) * \
                (self.D[self.hash_d[key_dict["t"]]] - \
                 (self.x[self.hash_x[key_dict["hj"]]] * self.y[self.hash_y[key_dict["hl"]]]) - \
                 (self.x[self.hash_x[key_dict["fk"]]] * self.y[self.hash_y[key_dict["fl"]]]) - \
                 (self.x[self.hash_x[key_dict["gi"]]] * self.y[self.hash_y[key_dict["gl"]]]))
    return numerator

def r_denominator(self,key_dict):
    denominator = self.tau[self.hash_t[key_dict["t"]]] * self.ss[self.hash_s[key_dict["s"]]]
    return denominator

def collect_all_keys(self,k_triple):
    #only get keys from tau and split them accordingly
    tau_keys = [x for x in self.hash_t.keys() if k_triple in x]
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