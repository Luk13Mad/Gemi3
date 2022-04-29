import pandas as pd
import numpy as np
import sys
import copy


def calc_observed_LFC(raw_data,ETP,replicates,normalize,constant=32):
    '''
    calculates the LFC score between early and late timepoint

    :param raw_data: df with raw count data, optionally normalized to million reads
    :param ETP: list with early timepoint column names
    :param replicates: replicate annotation
    :param normalize: bool, True or False
    :param constant: constant to add during LFC calculation, default 32
    :return: df with LFC table
    '''
    #raw_data="/Users/lukasmadenach/Desktop/Paper/Zhou2020/gemi3_on_zhou_data/count_table_2d_zhou.txt"
    #ETP=["Cell_D15.rep1","Cell_D15.rep2" ]
    #annotations="/Users/lukasmadenach/Desktop/Paper/Zhou2020/gemi3_on_zhou_data/annotation_2d_zhou.txt"
    #replicates="/Users/lukasmadenach/Desktop/Paper/Zhou2020/gemi3_on_zhou_data/replicates_2D_zhou.txt"
    raw_df=pd.read_csv(raw_data,sep="\t",index_col=0)
    rep_anno=pd.read_csv(replicates,sep="\t") #concatenated sample and rep must be separated by "."

    if normalize:
        total_counts=raw_df.sum(skipna=True).sum(skipna=True)
        scale=total_counts/raw_df.shape[1]
        raw_df=raw_df.apply(lambda x: ((x / x.sum())*scale)+constant, axis=0)

    raw_df=raw_df.apply(lambda x: np.log2(x)-np.log2(np.nanmedian(x)),axis=0)
    ETP_col=raw_df.loc[:,ETP].mean(axis=1)
    ETP_col.rename(ETP[0].split(".")[0],inplace=True)
    LFC_out=pd.DataFrame(ETP_col)
    rep_anno = rep_anno.loc[~rep_anno["colname"].isin(ETP), :]

    for s in set(rep_anno["samplename"]):
        mask=rep_anno.loc[rep_anno["samplename"]==s,:]["colname"]
        sample_col = raw_df.loc[:, list(mask)].mean(axis=1)
        sample_col.rename(s,inplace=True)
        LFC_out[s]=sample_col

    for s in LFC_out.columns:
        LFC_out[s]=LFC_out[s].subtract(ETP_col)
    LFC_out=LFC_out.drop(ETP[0].split(".")[0],axis=1)
    return LFC_out

class G3D():

    def __init__(self,LFC_table, annotation_table,samples,verbose=False):
        '''
        reads file with LFC table, also needs annotation table example files in Gemi3/data

        :param LFC_table: table wiht logfoldchange data, each column is one sample, index is guide RNA combinations, example in Gemi3/data
        :param annotation_table: table with guide annotation, used for matching guide RNA to genename, example in Gemi3/data
        :return: two values, LFC table values as pandas dataframe (df) and as tensorflow tf.connstant (observed)
        '''
        self.LFC_table=LFC_table
        self.samples=samples
        if "_" in samples:
            raise ValueError(" _ must not be in the sample names")
        self.annotation_table=annotation_table
        self.verbose=verbose
        if verbose: print("Loading LFC table from file.",sys.stdout)
        self.full_df = pd.read_csv(LFC_table, sep="\t",index_col=0)  # table with log fold changes between early and late timepoint
        if verbose: print("Loading annotation table from file.", sys.stdout)
        self.raw_anno_table,self.anno_dict_1, self.anno_dict_2, self.anno_dict_3 = self.read_annotation(annotation_table)
        if verbose: print("Merging annotation with LFC table from file.", sys.stdout)
        self.gene_pairs1 = [self.gRNApair_to_genenames(name, self.anno_dict_1, self.anno_dict_2, self.anno_dict_3) for name in self.full_df.index]
        self.full_df["first_gene"] = [x[0] for x in self.gene_pairs1]
        self.full_df["second_gene"] = [x[1] for x in self.gene_pairs1]
        self.full_df["third_gene"] = [x[2] for x in self.gene_pairs1]
        self.raw_data = self.full_df.loc[:, self.samples].values

    def __copy__(self):
        return G3D(self.LFC_table, self.annotation_table, self.samples,self.verbose)

    def __deepcopy__(self,memo):
        return G3D(copy.deepcopy(self.LFC_table,memo),copy.deepcopy(self.annotation_table,memo),
                   copy.deepcopy(self.samples,memo),copy.deepcopy(self.verbose,memo))

    def read_annotation(self,annotation_table):
        '''
        function for reading annotation table

        :param annotation_table: table with guide annotation, used for matching guide RNA to genename, example in Gemi3/data
        :return: three dicts with guide RNA as key and gene name as value, one dict for each spot in 3D vector
        '''
        annotation = pd.read_csv(annotation_table, sep="\t",index_col=0)
        if not all([x in annotation.columns for x in ["guide_1","guide_2","guide_3","gene_1","gene_2","gene_3"]]):
            print("There should be the columns guide_1,guide_2,guide_3,gene_1,gene_2 and gene_3",file=sys.stderr)
            raise ValueError
        _anno_dict_1 = {k: v for k, v in zip(annotation["guide_1"], annotation["gene_1"])}
        _anno_dict_2 = {k: v for k, v in zip(annotation["guide_2"], annotation["gene_2"])}
        _anno_dict_3 = {k: v for k, v in zip(annotation["guide_3"], annotation["gene_3"])}
        return annotation,_anno_dict_1, _anno_dict_2, _anno_dict_3

    def gRNApair_to_genenames(self,name, dict1, dict2, dict3):
        '''
        function to replace guide RNA with genenames

        :param name: concatenated guide RNAs, separated by ";", usually from index of dataframe
        :param dict1: dict with guide RNA as key and genenames as value
        :param dict2: dict with guide RNA as key and genenames as value
        :param dict3: dict with guide RNA as key and genenames as value
        :return: three values, input string split into separate guide RNA and guide RNA replaced by genename
        '''
        p1, p2, p3 = name.split(";")
        if p1 in dict1.keys() and p2 in dict2.keys() and p3 in dict3.keys():
            return dict1[p1], dict2[p2], dict3[p3]
        elif p2 in dict1.keys() and p1 in dict2.keys() and p3 in dict3.keys():
            return dict1[p2], dict2[p1], dict3[p3]
        elif p1 in dict1.keys() and p3 in dict2.keys() and p2 in dict3.keys():
            return dict1[p1], dict2[p3], dict3[p2]
        elif p3 in dict1.keys() and p2 in dict2.keys() and p1 in dict3.keys():
            return dict1[p3], dict2[p2], dict3[p1]


    def set_raw_df(self,new_data):
        '''
        function sets new_data in self.raw_data

        :param new_data: pandas dataframe with new values
        :return: -
        '''
        if isinstance(new_data, pd.DataFrame):
            if self.raw_data.shape == new_data.values.shape:
                self.raw_data = new_data.values
                self.set_full_df(self.raw_data)
            else:
                print("Dimensions don't match", sys.stderr)

    def set_full_df(self,new):
        '''
        function sets new in self.full_df together with columns first_gene, second_gene and third_gene

        :param new:  dataframe
        :return: -
        '''
        if self.full_df.loc[:,self.samples].values.shape == new.shape:
            self.full_df=pd.DataFrame(data=new,index=self.full_df.index,columns=self.samples)
            self.full_df["first_gene"] = [x[0] for x in self.gene_pairs1]
            self.full_df["second_gene"] = [x[1] for x in self.gene_pairs1]
            self.full_df["third_gene"] = [x[2] for x in self.gene_pairs1]
        else:
            print("Dimensions don't match",sys.stderr)