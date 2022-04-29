import numpy as np
import sys
from Gemi3.pyCAVI.util import find_s_keys,find_r_keys,gene_to_guides_fast

def y_loop(self):
    if self.verbose:
        print("Updating y values.",file=sys.stdout)
    for gl in self.hash_y.keys():
        indices=collect_all_keys(self,gl)
        numerator_sum=np.sum(np.asarray(list(map(lambda x: y_numerator(self,x),indices))))
        denominator_sum=np.sum(np.asarray(list(map(lambda x: y_denominator(self,x),indices))))
        new_ygl=((self.prio_tau_y*self.prio_mu_y)+numerator_sum)/(self.prio_tau_y+denominator_sum)
        if self.debug: print(f"Old ygl {self.y[self.hash_y[gl]]}.", sys.stdout)
        self.y[self.hash_y[gl]] = new_ygl
        if self.debug: print(f"New ygl {self.y[self.hash_y[gl]]}.", sys.stdout)


def yy_loop(self):
    if self.verbose:
        print("Updating yy values.",file=sys.stdout)
    for gl in self.hash_y.keys():
        indices=collect_all_keys(self,gl)
        numerator_sum=np.sum(np.asarray(list(map(lambda x: y_numerator(self,x),indices))))
        denominator_sum=np.sum(np.asarray(list(map(lambda x: y_denominator(self,x),indices))))
        new_yygl=((self.prio_tau_yy*self.prio_mu_yy)+numerator_sum)/(self.prio_tau_yy+denominator_sum)
        if self.debug: print(f"Old yygl {self.yy[self.hash_y[gl]]}.", sys.stdout)
        self.yy[self.hash_y[gl]] = new_yygl
        if self.debug: print(f"New yygl {self.yy[self.hash_y[gl]]}.", sys.stdout)

def y_numerator(self,key_dict):
    numerator=(self.tau[self.hash_t[key_dict["t"]]]*self.x[self.hash_x[key_dict["gi"]]])*\
    (self.D[self.hash_d[key_dict["t"]]]-\
    (self.x[self.hash_x[key_dict["hj"]]]*self.y[self.hash_y[key_dict["hl"]]])-\
    (self.x[self.hash_x[key_dict["fk"]]]*self.y[self.hash_y[key_dict["fl"]]])-\
    (self.r[self.hash_r[key_dict["r"]]]*self.s[self.hash_s[key_dict["s"]]]))
    return numerator

def y_denominator(self,key_dict):
    denominator = self.tau[self.hash_t[key_dict["t"]]] * self.xx[self.hash_x[key_dict["gi"]]]
    return denominator

def collect_all_keys(self,gl_extern):
    #only get keys from tau and split them accordingly
    gene,sample=gl_extern.split("_")
    gi_extern=gene_to_guides_fast(self,gene)
    tau_keys = [x for x in self.hash_t.keys() if any(list(map(lambda d: d in x,gi_extern)))]
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