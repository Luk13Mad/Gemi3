import numpy as np

def gene_to_guides(self,gene):
    sub = self.G3D.raw_anno_table.__deepcopy__()
    sub = sub.loc[(sub["gene_1"] == gene), :]
    sub = sub.loc[(sub["gene_2"] == gene), :]
    sub = sub.loc[(sub["gene_3"] == gene), :]
    return list(np.unique(sub[["guide_1", "guide_2", "guide_3"]]))

def gene_to_guides_fast(self,gene):
    if gene in self.gene_to_guide_dict:
        return self.gene_to_guide_dict[gene]
    else:
        return None

def guide_to_gene(self, guide):
    if guide in self.G3D.anno_dict_1:
        return self.G3D.anno_dict_1[guide]
    elif guide in self.G3D.anno_dict_2:
        return self.G3D.anno_dict_2[guide]
    elif guide in self.G3D.anno_dict_3:
        return self.G3D.anno_dict_3[guide]


def find_s_keys(self,gl,hl,fl,l): #for a given triple finds the s keys
    if "_".join([gl,hl,fl,l]) in self.hash_s:
        return "_".join([gl,hl,fl,l])
    else:
        for k in self.hash_s.keys():
            first,second,third,sample=k.split("_")
            if sample != l:
                continue
            elif gl == first and hl == second and fl == third:#123
                return k
            elif gl == second and hl == first and fl == third:#213
                return k
            elif gl == first and hl == third and fl == second:#132
                return k
            elif gl == second and hl == third and fl == first:#231
                return k
            elif gl == third and hl == first and fl == second:#312
                return k
            elif gl == third and hl == second and fl == first:#321
                return k
        return None

def find_r_keys(self,gi,hj,fk):#for a given triple finds the r keys
    if "_".join([gi,hj,fk]) in self.hash_r:
        return "_".join([gi,hj,fk])
    else:
        for k in self.hash_r.keys():
            first,second,third=k.split("_")
            if gi == first and hj == second and fk == third:#123
                return k
            elif gi == second and hj == first and fk == third:#213
                return k
            elif gi == first and hj == third and fk == second:#132
                return k
            elif gi == second and hj == third and fk == first:#231
                return k
            elif gi == third and hj == first and fk == second:#312
                return k
            elif gi == third and hj == second and fk == first:#321
                return k
        return None

def decreasing(self,L):
    '''
    Checks if list L is decreasing. NOT strictly decreasing.
    :param L:
    :return:
    '''
    return all(x>=y for x, y in zip(L, L[1:]))