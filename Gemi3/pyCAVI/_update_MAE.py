import numpy as np
from Gemi3.pyCAVI.util import find_s_keys,find_r_keys

def calc_MAE(self):
    indices=collect_all_keys(self)
    MAE=np.mean(np.abs(np.asarray(list(map(lambda x: MAE_one_triple(self,x),indices)))))
    return MAE


def MAE_one_triple(self,key_dict):
    out=self.D[self.hash_d[key_dict["t"]]]- \
        self.x[self.hash_x[key_dict["gi"]]] * self.y[self.hash_y[key_dict["gl"]]]- \
        self.x[self.hash_x[key_dict["hj"]]] * self.y[self.hash_y[key_dict["hl"]]]- \
        self.x[self.hash_x[key_dict["fk"]]] * self.y[self.hash_y[key_dict["fl"]]]- \
        self.r[self.hash_r[key_dict["r"]]]*self.s[self.hash_s[key_dict["s"]]]
    return out

def collect_all_keys(self):
    #only get keys from tau and split them accordingly
    tau_keys = [x for x in self.hash_t.keys()]
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