import sys
from Gemi3.pyCAVI.util import find_s_keys,find_r_keys

def tau_loop(self):
    if self.verbose:
        print("Updating tau values.",file=sys.stdout)
    for tau_triple in self.hash_t.keys():
        indices=collect_all_keys(self,tau_triple)
        new_alpha = tau_alpha(self,indices)
        self.alpha[self.hash_t[tau_triple]] = new_alpha
        new_beta=tau_beta(self,indices)
        self.beta[self.hash_t[tau_triple]] = new_beta
        new_tau=new_alpha/new_beta
        if self.debug: print(f"Old tau {self.tau[self.hash_t[tau_triple]]}.", sys.stdout)
        self.tau[self.hash_t[tau_triple]] = new_tau
        if self.debug: print(f"New s {self.tau[self.hash_t[tau_triple]]}.", sys.stdout)

def tau_alpha(self,key_dict):
    return self.alpha[self.hash_t[key_dict["t"]]]+0.5

def tau_beta(self,key_dict):
    beta=self.beta[self.hash_t[key_dict["t"]]]
    D2=0.5 * self.D[self.hash_d[key_dict["t"]]]**2
    Dxgiygl=self.D[self.hash_d[key_dict["t"]]]*self.x[self.hash_x[key_dict["gi"]]]*self.y[self.hash_y[key_dict["gl"]]]
    Dxhjyhl=self.D[self.hash_d[key_dict["t"]]]*self.x[self.hash_x[key_dict["hj"]]]*self.y[self.hash_y[key_dict["hl"]]]
    Dxfkyfl=self.D[self.hash_d[key_dict["t"]]]*self.x[self.hash_x[key_dict["fk"]]]*self.y[self.hash_y[key_dict["fl"]]]
    Drs=self.D[self.hash_d[key_dict["t"]]]*self.r[self.hash_r[key_dict["r"]]]*self.s[self.hash_s[key_dict["s"]]]
    xgi2ygl2=self.xx[self.hash_x[key_dict["gi"]]]*self.yy[self.hash_y[key_dict["gl"]]]
    xgiyglxhjyhl=self.x[self.hash_x[key_dict["gi"]]]*self.y[self.hash_y[key_dict["gl"]]]*self.x[self.hash_x[key_dict["hj"]]]*self.y[self.hash_y[key_dict["hl"]]]
    xgiyglxfkyfl=self.x[self.hash_x[key_dict["gi"]]]*self.y[self.hash_y[key_dict["gl"]]]*self.x[self.hash_x[key_dict["fk"]]]*self.y[self.hash_y[key_dict["fl"]]]
    xgiyglrs=2*self.x[self.hash_x[key_dict["gi"]]]*self.y[self.hash_y[key_dict["gl"]]]*self.r[self.hash_r[key_dict["r"]]]*self.s[self.hash_s[key_dict["s"]]]
    xhj2yhl2=0.5*self.xx[self.hash_x[key_dict["hj"]]]*self.yy[self.hash_y[key_dict["hl"]]]
    xhjyhlxfkyfl=self.x[self.hash_x[key_dict["hj"]]]*self.y[self.hash_y[key_dict["hl"]]]*self.x[self.hash_x[key_dict["fk"]]]*self.y[self.hash_y[key_dict["fl"]]]
    xhjyhlrs=self.x[self.hash_x[key_dict["hj"]]]*self.y[self.hash_y[key_dict["hl"]]]*self.r[self.hash_r[key_dict["r"]]]*self.s[self.hash_s[key_dict["s"]]]
    r2s2=0.5*self.rr[self.hash_r[key_dict["r"]]]*self.ss[self.hash_s[key_dict["s"]]]
    out=beta+D2-Dxgiygl-Dxhjyhl-Dxfkyfl-Drs+xgi2ygl2+xgiyglxhjyhl+xgiyglxfkyfl+xgiyglrs+xhj2yhl2+xhjyhlxfkyfl+xhjyhlrs+r2s2
    return out


def collect_all_keys(self,tau_triple):
    gi,hj,fk,l = tau_triple.split("_")
    gl = self.guide_to_gene(gi)
    hl = self.guide_to_gene(hj)
    fl = self.guide_to_gene(fk)
    s = find_s_keys(self,gl,hl,fl,l)
    r = find_r_keys(self,gi,hj,fk)
    out={"t":tau_triple,"gi":gi,"hj":hj,"fk":fk,"l":l,
                    "gl":"_".join([gl,l]),
                    "hl":"_".join([hl,l]),
                    "fl":"_".join([hl,l]),
                    "r":r,
                    "s":s}
    return out