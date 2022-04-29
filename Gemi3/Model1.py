import numpy as np
import sys
from . import g3datahandling as g3d
import itertools as itt
import math
from Gemi3.rsCAVI import start_inference


class Model1:
    from Gemi3.pyCAVI._update_x import x_loop,xx_loop #modules with coordinate update equations
    from Gemi3.pyCAVI._update_y import y_loop,yy_loop
    from Gemi3.pyCAVI._update_r import r_loop,rr_loop
    from Gemi3.pyCAVI._update_s import s_loop,ss_loop
    from Gemi3.pyCAVI._update_tau import tau_loop
    from Gemi3.pyCAVI._update_MAE import calc_MAE
    from Gemi3.pyCAVI.util import gene_to_guides,decreasing,guide_to_gene,gene_to_guides_fast  #guide_to_gene and gene_to_guides_fast functions needed as method for other scripts to access

    def __init__(self,G3D,ncgenes,verbose=False,debug=False,force=False,MAX_ITERATIONS=30,MAE_THRESHOLD=0.001,
                 prio_mu_x=1.0,prio_tau_x=1.0,prio_mu_xx=1.0,prio_tau_xx=1.0,
                 prio_mu_y=0.0,prio_tau_y=10.0,prio_mu_yy=0.0,prio_tau_yy=10.0,
                 prio_mu_r=1.0,prio_tau_r=1.0,prio_mu_rr=1.0,prio_tau_rr=1.0,
                 prio_mu_s=0.0,prio_tau_s=10.0,prio_mu_ss=0.0,prio_tau_ss=10.0):
        if not isinstance(G3D,g3d.G3D): # check if data is in correct data handling object
            raise TypeError("G3D object needed as input")
        self.G3D = G3D #set field with full datahandling object
        self.debug=debug #set bool for debug mode
        if self.debug: #if debug mode verbose is automatically true
            print("!!!!MODEL initialized in DEBUG mode!!!!\n!!!!Very verbose output will be generated!!!!")
            self.verbose=True
        else:
            self.verbose = verbose
        self.force=force #set bool for force, if force==true CAVI continues even if MAE increases after 1 step
        self.full_df = G3D.full_df #set field with pandas dataframe holding LFC and gene annotation
        self.MAX_ITERATIONS=MAX_ITERATIONS #set field with maximum number of iterations for CAVI
        self.MAE_THRESHOLD=MAE_THRESHOLD #set field with threshold. if change in MAE smaller than threshold CAVI stops
        self.samples = G3D.samples #set field with list of sample names
        self.is_primed = False #set state
        self.is_inferred = False #set inference state
        self.D = G3D.raw_data #set numpy array with LFC values
        if isinstance(ncgenes,list): #check if nc-genes were supplied and are correct type
            if len(ncgenes)==0:
                raise TypeError("ncgenes must be supplied")
            self.ncgenes=ncgenes
        else:
            raise TypeError("Ncgenes must be supplied as list of strings.")
        self.prio_mu_x=prio_mu_x #x and xx prio
        self.prio_tau_x=prio_tau_x
        self.prio_mu_xx = prio_mu_xx
        self.prio_tau_xx = prio_tau_xx
        self.prio_mu_y=prio_mu_y # y and yy prio
        self.prio_tau_y=prio_tau_y
        self.prio_mu_yy = prio_mu_yy
        self.prio_tau_yy = prio_tau_yy
        self.prio_mu_r=prio_mu_r # r and rr prio
        self.prio_tau_r=prio_tau_r
        self.prio_mu_rr = prio_mu_rr
        self.prio_tau_rr = prio_tau_rr
        self.prio_mu_s=prio_mu_s # s and ss prio
        self.prio_tau_s=prio_tau_s
        self.prio_mu_ss = prio_mu_ss
        self.prio_tau_ss = prio_tau_ss
        self.ncguides = list(itt.chain(*[self.gene_to_guides(x) for x in self.ncgenes])) #convert supplied ncgenes to guides and set
        if len(self.ncguides)==0:
            raise TypeError("Specified ncgenes not present")
        ncgenes_mask = self.full_df.first_gene.isin(self.ncgenes) | self.full_df.second_gene.isin(
            self.ncgenes) | self.full_df.third_gene.isin(self.ncgenes) #bool mask for removing rows containg ncgenes from df
        self.ncgenes_df = self.full_df.loc[ncgenes_mask, :] #set pandas df with only ncgenes
        self.full_df = self.full_df.loc[~ncgenes_mask, :] #remove all rows containing at least 1 nc gene
        only_two_target_genes_mask = self.full_df.apply(lambda x: x["first_gene"]==x["second_gene"] or x["first_gene"]==x["third_gene"] or x["second_gene"]==x["third_gene"],axis=1) #set mask to filter rows were only two different genes are targeted
        if self.verbose:
            print("Removing rows where only two different genes are targeted",file=sys.stdout)
            print(f"Removing {np.sum(only_two_target_genes_mask)} rows.",file=sys.stdout)
        self.full_df = self.full_df.loc[~only_two_target_genes_mask,:]
        if self.full_df.empty:
            raise ValueError("All rows only contain two different genes.")
        self.hash_x=dict()#initialze hash dict empty, later used in prime() method
        self.hash_y=dict()
        self.hash_r=dict()
        self.hash_s=dict()
        self.hash_t=dict()
        self.hash_d=dict()
        self.MAE=[] #set MAE as empty list
        unique_guides_tmp = list(map(lambda x: x.split(";"), self.full_df.index))
        self.unique_guides = [x for x in set(itt.chain(*unique_guides_tmp)) if x not in self.ncguides] #set list with unique guides removing ncguides
        del unique_guides_tmp
        self.unique_genes = [x.replace(" ","") for x in np.unique(self.full_df.loc[:, ["first_gene", "second_gene", "third_gene"]].values) if x not in self.ncgenes] #set list with unique genes (not including ncgenes) and removing whitespace from gene names
        self.gene_to_guide_dict={k:self.gene_to_guides(k) for k in self.unique_genes}
        self.guide_to_gene_dict=dict() #set guide to gene dict
        self.guide_to_gene_dict.update(G3D.anno_dict_1)
        self.guide_to_gene_dict.update(G3D.anno_dict_2)
        self.guide_to_gene_dict.update(G3D.anno_dict_3)


    def prime(self):
        '''
        functions sets values from input data, build hash tables, initializes variables from control genes

        :return: -
        '''
        unique_guides_len = len(self.unique_guides)
        unique_genes_len = len(self.unique_genes)
        unique_threeway_combinations_guides = len(self.full_df.loc[:, ["first_gene", "second_gene", "third_gene"]].index)
        unique_threeway_combinations_genes = math.comb((unique_genes_len),(3))

        #D values
        self.D=self.G3D.raw_data
        self.build_hashes("D") #hashes for D same as for tau

        #x values
        self.x = np.random.normal(self.prio_mu_x,np.reciprocal(self.prio_tau_x),(unique_guides_len))
        self.build_hashes("x")  # build hashes for x
        self.init_x_from_ncgenes() #init x values form ncgenes
        self.xx = self.x**2

        #y values
        self.y = np.random.normal(self.prio_mu_y,np.reciprocal(self.prio_tau_y),(unique_genes_len,len(self.samples)))
        self.build_hashes("y")  # build hashes for y
        self.init_y_from_ncgenes() #init y values form ncgenes
        self.yy = self.y ** 2

        #r values
        self.r = np.random.normal(self.prio_mu_r,np.reciprocal(self.prio_tau_r),(unique_threeway_combinations_guides))
        self.build_hashes("r")  # build hashes for r
        self.rr = self.r**2

        #s values
        self.s = np.random.normal(self.prio_mu_s,np.reciprocal(self.prio_tau_s),(unique_threeway_combinations_genes,len(self.samples)))
        self.build_hashes("s")  # build hashes for s
        self.ss = self.s**2

        #tau values
        self.tau = np.random.gamma(0.5,1.0,(len(self.full_df.index),len(self.samples)))
        self.build_hashes("t")  # build hashes for tau
        self.alpha = np.full((len(self.full_df.index),len(self.samples)),0.5)
        self.beta = np.full((len(self.full_df.index),len(self.samples)),1.0) #TODO write init function for beta

        #MAE values
        self.MAE.append(self.calc_MAE())
        if self.verbose:
            print(f"Starting MAE at {self.MAE[-1]}")
        self.is_primed=True

    def build_hashes(self,variable):
        '''
        Builds and stores hashes needed for faster CAVI

        :param variable:
        :return:
        '''
        unique_threeway_combinations_guides = self.full_df.loc[:, ["first_gene", "second_gene", "third_gene"]].index
        unique_threeway_combinations_genes = math.comb(len(self.unique_genes), (3))

        if variable=="x": #initialization removed
            if self.verbose: print("Building hashtable for x and xx",file=sys.stdout)
            if len(self.unique_guides) != len(self.x):
                #print("Error initializing x and xx",file=sys.stderr)
                raise ValueError("Error initializing x and xx")
            for guide,idx in zip(self.unique_guides,range(len(self.x))):
                self.hash_x.update({guide:idx})
            if self.debug: print(self.hash_x.keys(),sys.stdout)

        elif variable == "y": #initialization removed
            if self.verbose: print("Building hashtable for y and yy",file=sys.stdout)
            if self.y.shape != (len(self.unique_genes),len(self.samples)):
                #print("Error initializing y and yy",file=sys.stderr)
                raise ValueError("Error initializing y and yy")
            num_rows, num_cols = self.y.shape
            for guide,idx1 in zip(self.unique_genes,range(num_rows)):
                for sample,idx2 in zip(self.samples,range(num_cols)):
                    self.hash_y.update({"_".join([guide,sample]):(idx1,idx2)})
            if self.debug: print(self.hash_y.keys(), sys.stdout)

        elif variable == "r":#all guide combinations in specified order, initialization genes are removed
            if self.verbose: print("Building hashtable for r and rr", file=sys.stdout)
            if len(self.r) != len(unique_threeway_combinations_guides):
                #print("Error initializing r and rr", file=sys.stderr)
                raise ValueError("Error initializing r and rr")
            sorted_threeway_guides=[x.replace(";","_") for x in unique_threeway_combinations_guides]
            for guide,idx in zip(sorted_threeway_guides,range(len(self.r))):
                self.hash_r.update({guide:idx})
            if self.debug: print(self.hash_r.keys(), sys.stdout)

        elif variable == "s": #should not care about order of genes in combination, initialization genes are removed
            if self.verbose: print("Building hashtable for s and ss",file=sys.stdout)
            if self.s.shape != (unique_threeway_combinations_genes,len(self.samples)):
                #print("Error initializing s and ss", file=sys.stderr)
                raise ValueError("Error initializing s and ss")
            sorted_threeway_genes = np.sort(list(self.unique_genes))
            sorted_threeway_genes = itt.combinations(sorted_threeway_genes,3)
            num_rows, num_cols = self.s.shape
            for gene,idx1 in zip(sorted_threeway_genes,range(num_rows)):
                for sample,idx2 in zip(self.samples,range(num_cols)):
                    self.hash_s.update({"_".join(["_".join(gene),sample]):(idx1,idx2)})
            if self.debug: print(self.hash_r.keys(), sys.stdout)

        elif variable == "t":
            if self.verbose: print("Building hashtable for tau",file=sys.stdout)
            if self.tau.shape != (len(self.full_df.index),len(self.samples)):
                #print("Error initializing tau", file=sys.stderr)
                raise ValueError("Error initializing tau")
            num_rows, num_cols = self.tau.shape
            sorted_threeway_guides = [x.replace(";","_") for x in self.full_df.index]
            for guide,idx1 in zip(sorted_threeway_guides,range(num_rows)):
                for sample,idx2 in zip(self.samples,range(num_cols)):
                    self.hash_t.update({"_".join([guide,sample]):(idx1,idx2)})
            if self.debug: print(self.hash_t.keys(), sys.stdout)

        elif variable == "D":
            if self.verbose: print("Building hashtable for D",file=sys.stdout)
            num_rows, num_cols = self.D.shape
            sorted_threeway_guides = [x.replace(";","_") for x in self.full_df.index]
            for guide,idx1 in zip(sorted_threeway_guides,range(num_rows)):
                for sample,idx2 in zip(self.samples,range(num_cols)):
                    self.hash_d.update({"_".join([guide,sample]):(idx1,idx2)})
            if self.debug: print(self.hash_t.keys(), sys.stdout)

        else:
            print("Build hashes failed, should not happen, check source code",file=sys.stderr)
            raise KeyError("Wrong key for building hashes")

    def start_inference_python(self): #inferenc method, each variable gets updated while others are constant
        '''
        functions runs CAVI algorithm implemented with python
        !!!ATTENTION!!! faults are in this implementation of CAVI, use rsCAVI

        :return:
        '''
        print("!!!ATTENTION!!!\n!!!DEPRECATED!!!\nuse rsCAVI instead of this function",file=sys.stderr)
        if not self.is_primed:
            raise TypeError("Model must be primed first.")
        counter=0
        while counter <= self.MAX_ITERATIONS:
            self.x_loop()
            self.xx_loop()
            self.y_loop()
            self.yy_loop()
            self.r_loop()
            self.rr_loop()
            self.s_loop()
            self.ss_loop()
            self.tau_loop()
            self.MAE.append(self.calc_MAE())
            counter +=1
            if self.verbose:
                print(f"Iteration {counter} complete.",file=sys.stdout)
                print(f"MAE at {self.MAE[-1]}.",file=sys.stdout)
            if not self.force:
                if self.MAE[-1]>self.MAE[-2]:
                    print("Breaking CAVI: MAE increased with step. \n Init model with force=True if you want to keep going.", sys.stdout)
                    break
            if np.abs(self.MAE[-1]-self.MAE[-2])<self.MAE_THRESHOLD:
                if self.verbose: print("Breaking CAVI: MAE_THRESHOLD reached.",sys.stdout)
                break
            if self.MAE[-1] < 0:
                print("Breaking CAVI: MAE turned negative.", sys.stdout)
                break

    def start_inference_rust(self):
        '''
        functions runs CAVI algorithm implemented with rust

        :return:
        '''
        if not self.is_primed:
            raise TypeError("Model must be primed first.")
        if self.verbose:
            print("Handoff to Rust.", sys.stdout)
        result_from_rust = start_inference(self.verbose, self.debug, self.force, self.MAX_ITERATIONS, self.MAE_THRESHOLD,
                                   self.ncgenes, self.ncguides, self.samples, self.MAE,
                                   self.unique_guides, self.unique_genes, self.guide_to_gene_dict,
                                   self.gene_to_guide_dict,
                                   self.prio_mu_x, self.prio_tau_x, self.prio_mu_xx, self.prio_tau_xx,
                                   self.prio_mu_y, self.prio_tau_y, self.prio_mu_yy, self.prio_tau_yy,
                                   self.prio_mu_r, self.prio_tau_r, self.prio_mu_rr, self.prio_tau_rr,
                                   self.prio_mu_s, self.prio_tau_s, self.prio_mu_ss, self.prio_tau_ss,
                                   self.hash_x, self.x.tolist(), self.xx.tolist(), self.x.shape[0],
                                   self.hash_y, self.y.flatten().tolist(), self.yy.flatten().tolist(),
                                   self.y.shape[0], self.y.shape[1],
                                   self.hash_r, self.r.tolist(), self.rr.tolist(), self.r.shape[0],
                                   self.hash_s, self.s.flatten().tolist(), self.ss.flatten().tolist(),
                                   self.s.shape[0], self.s.shape[1],
                                   self.hash_t, self.tau.flatten().tolist(), self.tau.shape[0],
                                   self.tau.shape[1],
                                   self.alpha.flatten().tolist(), self.beta.flatten().tolist(),
                                   self.full_df[self.samples].to_numpy().flatten().tolist())
        if self.verbose:
            print("Result returned from Rust. Checking...", sys.stdout)
        if not self.decreasing(result_from_rust[0]):
            print("MAE increased with step. \n Discarding results.", sys.stdout)
            return
        if self.verbose:
            print("Check completed. Setting values.", sys.stdout)
        nrow_y = self.y.shape[0]
        ncol_y = self.y.shape[1]
        nrow_s = self.s.shape[0]
        ncol_s = self.s.shape[1]
        nrow_tau = self.tau.shape[0]
        ncol_tau = self.tau.shape[1]

        self.MAE = result_from_rust[0]
        self.x = np.array(result_from_rust[1])
        self.xx = np.array(result_from_rust[2])
        self.y = np.array(result_from_rust[3]).reshape((nrow_y,ncol_y))
        self.yy = np.array(result_from_rust[4]).reshape((nrow_y,ncol_y))
        self.r = np.array(result_from_rust[5])
        self.rr = np.array(result_from_rust[6])
        self.s = np.array(result_from_rust[7]).reshape((nrow_s,ncol_s))
        self.ss = np.array(result_from_rust[8]).reshape((nrow_s,ncol_s))
        self.tau = np.array(result_from_rust[9]).reshape((nrow_tau,ncol_tau))
        self.alpha = np.array(result_from_rust[10]).reshape((nrow_tau,ncol_tau))
        self.beta = np.array(result_from_rust[11]).reshape((nrow_tau,ncol_tau))
        self.is_inferred = True #set new inference state


    def init_x_from_ncgenes(self):
        '''
        calc initial values for x variables from specified control genes

        :return:
        '''
        counter = 0
        for guide in self.hash_x.keys():
            counter += 1
            mask = [x for x in self.ncgenes_df.index if guide in x]
            subdf = self.ncgenes_df.loc[mask, :] #subdf with guide in every row at least once
            two_nc_genes_mask = subdf.apply(
                lambda x: x["first_gene"] in self.ncgenes and x["second_gene"] in self.ncgenes or x[
                    "first_gene"] in self.ncgenes and x["third_gene"] in self.ncgenes or x[
                              "second_gene"] in self.ncgenes and x["third_gene"] in self.ncgenes, axis=1)
            subdf = subdf.loc[two_nc_genes_mask, :]
            if subdf.empty:
                continue
            init_x_value = np.median(subdf.loc[:,self.samples].values)
            self.x[self.hash_x[guide]] = init_x_value
        if counter != len(self.hash_x.keys()):
            print("WARNING: Not all genes paired with 2 ncgenes",sys.stderr)


    def init_y_from_ncgenes(self):
        '''
        calc initial values for y variables from specified control genes

        :return:
        '''
        counter = 0
        for gene_sample in self.hash_y.keys():
            counter += 1
            gene,sample=gene_sample.split("_")
            mask = self.ncgenes_df.first_gene.isin([gene]) | self.ncgenes_df.second_gene.isin([gene]) | self.ncgenes_df.third_gene.isin([gene])
            subdf = self.ncgenes_df.loc[mask,:] #subdf with gene at least once in every row
            two_nc_genes_mask = subdf.apply(
                lambda x: x["first_gene"] in self.ncgenes and x["second_gene"] in self.ncgenes or x["first_gene"] in self.ncgenes and x["third_gene"] in self.ncgenes or x[
                    "second_gene"] in self.ncgenes and x["third_gene"] in self.ncgenes, axis=1)
            subdf =subdf.loc[two_nc_genes_mask,:]
            if subdf.empty:
                continue
            init_y_value=np.median(subdf.loc[:,sample].values)
            self.y[self.hash_y[gene_sample]]=init_y_value
        if counter != len(self.hash_y.keys()):
            print("WARNING: Not all genes paired with 2 ncgenes",sys.stderr)

