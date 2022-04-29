import numpy as np
import sys
from Gemi3.pyCAVI.util import find_s_keys,find_r_keys

def x_loop(self):
    if self.verbose:
        print("Updating x values",file=sys.stdout)
    for gi in self.unique_guides:
        indices=collect_all_keys(self,gi)
        numerator_sum=np.sum(np.asarray(list(map(lambda x: x_numerator(self,x),indices))))
        denominator_sum=np.sum(np.asarray(list(map(lambda x: x_denominator(self,x),indices))))
        new_xgi=((self.prio_tau_x*self.prio_mu_x)+numerator_sum)/(self.prio_tau_x+denominator_sum)
        if self.debug:print(f"Old xgi {self.x[self.hash_x[gi]]}.",sys.stdout)
        self.x[self.hash_x[gi]]= new_xgi
        if self.debug: print(f"New xgi {self.x[self.hash_x[gi]]}.", sys.stdout)

def xx_loop(self):
    if self.verbose:
        print("Updating xx values",file=sys.stdout)
    for gi in self.unique_guides:
        indices = collect_all_keys(self, gi)
        numerator_sum = np.sum(np.asarray(list(map(lambda x: x_numerator(self, x), indices))))
        denominator_sum = np.sum(np.asarray(list(map(lambda x: x_denominator(self, x), indices))))
        new_xxgi = ((self.prio_tau_xx * self.prio_mu_xx) + numerator_sum) / (self.prio_tau_xx + denominator_sum)
        if self.debug: print(f"Old xxgi {self.xx[self.hash_x[gi]]}.", sys.stdout)
        self.xx[self.hash_x[gi]] = new_xxgi
        if self.debug: print(f"New xxgi {self.xx[self.hash_x[gi]]}.", sys.stdout)

def x_numerator(self,key_dict):
    numerator=(self.tau[self.hash_t[key_dict["t"]]]*self.y[self.hash_y[key_dict["gl"]]])*\
    (self.D[self.hash_d[key_dict["t"]]]-\
    (self.x[self.hash_x[key_dict["hj"]]]*self.y[self.hash_y[key_dict["hl"]]])-\
    (self.x[self.hash_x[key_dict["fk"]]]*self.y[self.hash_y[key_dict["fl"]]])-\
    (self.r[self.hash_r[key_dict["r"]]]*self.s[self.hash_s[key_dict["s"]]]))
    return numerator

def x_denominator(self,key_dict):
    denominator = self.tau[self.hash_t[key_dict["t"]]] * self.yy[self.hash_y[key_dict["gl"]]]
    return denominator


def collect_all_keys(self,gi_extern):
    #only get keys from tau and split them accordingly
    tau_keys = [x for x in self.hash_t.keys() if gi_extern in x] #TODO check in "if x.startswith(gi)" is needed
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

