# import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import statsmodels.stats.multitest as multi
import os
import platform
if platform.system() == 'Windows':
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"                            # avoid warning message when using sklearn on Windows


class Score(object):
    """
    This class comprises functions for scoring triple gene interactions (strong interactions, lethality and recovery).
    """

    def __init__(self, model, nc_genes, pc_genes=None, pc_weight=0.5, verbose=False, output_dir=os.getcwd(),plotting=False):
        '''
        The inferred values from the class Model1 as well as positive and negative control genes are used as input.

        :param model: sample-independent effect (x) of guide i to gene g, sample-dependent individual effect (y),
         sample-dependent combinatorial effect (s), mean absolute error (MAE), and other parameters
        :param pc_genes: string or list of strings of positive control genes for removing interactions which
        involve individually lethal genes
        :param nc_genes: string or list of strings of negative control genes for building a Gaussian model
        :param pc_weight: weight of positive control gene
        :param verbose : sets verbosity, defaults to False
        :param output_dir : sets output directory, defaults to current working directory
        :param plotting : bool if plots should be generated, defaults to False
        :return: 4 values (three pandas dataframes with scores of triple gene interactions, cell line info)
        '''
        self.model = model                                                      # initialise output of class Model1
        self.verbose = verbose                                                  # initialise verbose for logging output
        self.plotting = plotting
        if os.path.exists(output_dir):                                          # check if supplied output dir exists
            self.output_dir = output_dir
        else:
            try:
                os.makedirs(output_dir) #try creating supplied output dir
                self.output_dir = output_dir
            except:
                print(f"Setting output directory to {os.getcwd()}")  # or use cwd as default
                self.output_dir = os.getcwd()
        if not self.model.is_inferred:
            raise TypeError("Run start_inference_rust() before calculating scores.")    # abort __init__
        elif self.model.is_inferred and self.verbose:                           # explicit opt-in for verbose output
            print("Model1 was successfully updated.")
        if pc_genes is None:
            self.pc_genes = []
        else:
            try:
                self.pc_genes =  self.check_input(pc_genes,"Positive control")   # check positive control genes
            except:
                raise
        if len(self.pc_genes) == 0:                                             # check if any pc genes were supplied
            print("No sensitive scores will be calculated, as no positive control gene was provided.")
        try:
            self.nc_genes = self.check_input(nc_genes, "Negative control")          # check negative control genes
        except:
            raise
        if pc_weight != 0.5:
            self.pc_weight = pc_weight                                          # initialise weight of pc gene
        else:
            self.pc_weight = pc_weight
            if self.verbose:
                print("Weight for positive control gene was set to 0.5 by default.")
        self.score_result = self.scoring()         # define output of class Score
        if self.plotting:
            self.plots = self.plotting_function()                                            # initialise drawing of plots


    def check_input(self, gene, name : str):
        #TODO are control genes in data?
        '''
        Function to check if positive and negative control genes were supplied.
        :param gene: name of variable of positive or negative control genes
        :param name: purpose of variable ("Positive control" or "Negative control")
        :return: If true: initialized variable, if false: error message
        '''
        if len(gene) == 0:                                                      # check if any nc genes were supplied
            raise TypeError(name + " genes must be supplied")
        elif isinstance(gene, str):                                        # type string and list are allowed
            self.gene = gene
        elif isinstance(gene, list):
            if all([isinstance(x,str) for x in gene]):
                self.gene = gene
        else:
            raise TypeError(name + " genes must be supplied as str or list of strings")
        return self.gene


    def new_dic(self, res, hashes):
        '''
        Function to modify the data structure.
        :param res: gene names
        :param hashes: inferred values
        :return: dictionary containing gene names as keys and inferred values as values
        '''
        #gene_list = [] if gene_list is None else gene_list
        #num_list = [] if num_list is None else num_list
        #val_list = [] if val_list is None else val_list
        dic = dict()

        for key, value in hashes.items():                                       # iterate through dictionary
            #gene_list.append(key)                                               # add genes as keys to list
            #num_list.append(value)                                              # add numbers of genes to list
            dic.update({key:res[value]})

        #for i in num_list:                                                      # i=0 or i=(0,0)
        #    val_list.append(res[i])                                             # add value of genes/guides to list

        #dic = dict(zip(gene_list, val_list))                                    # convert two np arrays to dict
        return dic


    def df_triple(self, gene_triple):
        '''
        Function to create a data frame where triple genes are mapped to their single genes.
        :param gene_triple: dictionary containing names of three genes as keys and inferred values as values
        :return: pandas dataframe with names of three genes as row names and three columns with single gene names
        '''
        df = pd.DataFrame.from_dict(gene_triple, orient='index')                # create df from dict using keys as rows
        df[["gene1", "gene2", "gene3"]] = ""                                    # add 3 new columns with empty str
        df = df.drop(0, axis=1)                                                 # delete column with s-values

        for rowname, col in df.iterrows():                                      # loop through rows and columns of df
            col[0] = rowname.split('_')[0]                                      # split row name by underscore
            col[1] = rowname.split('_')[1]                                      # and map 1st element to 1st col etc.
            col[2] = rowname.split('_')[2]

        # remove cell line info in row name & keep only unique rows (be independent from cell lines)
        df = df.rename(index=lambda x: x.rsplit('_', 1)[0]).drop_duplicates()
        return df


    def cell_info(self, y):
        '''
        Function to extract the information how many and which cell lines exist in data.
        :param y: dictionary containing name of one gene and cell line as keys and inferred values as values
        :return: list with names of cell lines as strings
        '''
        cell_list = []
        for key in y.keys():                                                               # iterate through names of keys
            cell = key.rsplit('_', 1)[1]                                            # keep only cell line info
            cell_list.append(cell)
        cell_list = list(set(cell_list))                                  # remove duplicate entries
        return cell_list


    def get_y(self, genes, y, cell_list):
        '''
        Function to create a numpy array containing the y values for calculating scores
        :param genes: names of genes at position 1, 2, or 3 (g,h,f)
        :param y: dictionary containing name of one gene and cell line as keys and inferred values as values
        :param cell_list: list with names of cell lines as strings
        :return: numpy array containing the y values of genes at position 1, 2, or 3 for calculating scores,
                 if multiple cell lines, the output is a stacked np array with one column per cell line
        '''
        lst = []

        if len(cell_list) == 1:
            for gene in genes:                                                  # iterate through names in array
                for key in sorted(y.keys()):                                                   # iterate through names of keys
                    key_n = key.rsplit('_', 1)[0]                               # remove info about cell line
                    if gene == key_n:                                           # if gene name matches key:
                        lst.append(y[key])                                      # - append list by value of key
            return np.expand_dims(np.array(lst), axis=1)                        # return array with y values

        elif len(cell_list) > 1:
            array = np.zeros((len(genes), 1))
            for i in cell_list:                                                 # iterate through cell lines
                y_cell = {k: v for k, v in y.items() if i in k}   # get dictionary with y values for cell line
                for gene in genes:
                    for key in sorted(y_cell.keys()):                                          # iterate through keys of cell line
                        key_n = key.rsplit('_', 1)[0]                               # remove info about cell line
                        if gene == key_n:                                       # if gene name matches key:
                            lst.append(y_cell[key])                             # - append list by value of key
                gene_values = np.expand_dims(np.array(lst), axis=1)          # add one dimension to reshape values array
                array = np.column_stack((array, gene_values))
                lst.clear()                                                     # empty list before starting loop again
            array = np.delete(array, 0, 1)                                      # delete first column containing zeros
            return array

        else:
            raise TypeError("Cell line information missing")                    # len(y) == 0

    def get_s(self,genes,s,cell_list):
        '''
        Function to create a numpy array containing the s values for calculating scores
        :param genes: names of genes in triple
        :param y: dictionary containing name of one triple and cell line as keys and inferred values as values
        :param cell_list: list with names of cell lines as strings
        :return: numpy array containing the s values of genes  for calculating scores,
                 if multiple cell lines, the output is a stacked np array with one column per cell line
        '''
        lst = []

        if len(cell_list) == 1:
            for gene in genes:  # iterate through names in array
                for key in sorted(s.keys()):  # iterate through names of keys
                    key_n = key.rsplit('_', 1)[0]  # remove info about cell line
                    if gene == key_n:  # if gene name matches key:
                        lst.append(s[key])  # - append list by value of key
            return np.expand_dims(np.array(lst), axis=1)  # return array with y values

        elif len(cell_list) > 1:
            array = np.zeros((len(genes), 1))
            for i in cell_list:  # iterate through cell lines
                s_cell = {k: v for k, v in s.items() if i in k}  # get dictionary with y values for cell line
                for gene in genes:
                    for key in sorted(s_cell.keys()):  # iterate through keys of cell line
                        key_n = key.rsplit('_', 1)[0]  # remove info about cell line
                        if gene == key_n:  # if gene name matches key:
                            lst.append(s_cell[key])  # - append list by value of key
                gene_values = np.expand_dims(np.array(lst), axis=1)  # add one dimension to reshape values array
                array = np.column_stack((array, gene_values))
                lst.clear()  # empty list before starting loop again
            array = np.delete(array, 0, 1)  # delete first column containing zeros
            return array

        else:
            raise TypeError("Cell line information missing")

    def calculate_scores(self, yg, yh, yf, sghf, mod_y):
        '''
        Function to compute the scores strong, sensitive lethality and sensitive recovery.
        :param yg: numpy array containing the y values of genes at position 1 (for every cell line)
        :param yh: numpy array containing the y values of genes at position 2 (for every cell line)
        :param yf: numpy array containing the y values of genes at position 3 (for every cell line)
        :param sghf: inferred values for triple gene interactions
        :param mod_y: dictionary containing name of one gene and cell line as keys and inferred values as values
        :return: numpy array with scores in correct order
        '''
        strong = np.abs(sghf) - np.abs(yg) - np.abs(yh) - np.abs(yf)            # Calculate strong score

        # Sensitive scores are only calculated if positive control gene was provided
        if len(self.pc_genes) > 0:
            sensitive_lethality_a = np.amin([yg, yh, yf],axis=0) - yg - yh - yf - sghf
            sensitive_recovery_a = yg + yh + yf + sghf - np.amin([yg, yh, yf],axis=0)

            # Update sensitive scores because one gene can either have a lethality or recovery effect
            pc_gene_keylist = [x for x in mod_y if x.split('_')[0] in self.pc_genes]
            pc_values = [mod_y[i] for i in pc_gene_keylist]
            pc_gene_mean = np.asarray(pc_values).mean()                         # compute mean of cell lines for pc gene
            threshold = pc_gene_mean * self.pc_weight

            if threshold >= 0:
                # Lethality: keep pairs if yg & yh & yf > threshold (if true = 1, if false = 0)
                lethality = np.where((yg >= threshold) & (yh >= threshold) & (yf >= threshold), 1, 0)
                # Recovery: keep pairs if one gene is smaller than threshold
                recovery = np.where((yg < threshold) | (yh < threshold) | (yf < threshold), 1, 0)
            else:
                lethality = np.where((yg < threshold) & (yh < threshold) & (yf < threshold), 1, 0)
                recovery = np.where((yg >= threshold) | (yh >= threshold) | (yf >= threshold), 1, 0)

            # Filter scores with mask lethality / recovery
            sensitive_lethality = np.where((lethality == 1), sensitive_lethality_a, np.nan)
            sensitive_recovery = np.where((recovery == 1), sensitive_recovery_a, np.nan)

        # Create sensitive scores with same shape as strong score but filled with zeros
        else:
            sensitive_lethality = np.zeros(np.shape(strong))
            sensitive_recovery = np.zeros(np.shape(strong))

        return strong, sensitive_lethality, sensitive_recovery


    def nc_values(self, keys_tgenes):
        '''
        Function to define scores including negative control genes
        :param keys_tgenes: Numpy array with names of all triple gene combinations
        :return: np array containing names of triple gene combinations which include at least one nc gene
        '''
        if isinstance(self.nc_genes, type('')):                                             # check if type str
            ncgenes_values = np.array([x for x in keys_tgenes if self.nc_genes in x])
        elif isinstance(self.nc_genes, type([])):                                           # check if type list
            ncgenes_all = []
            for index, gene in enumerate(self.nc_genes):
                ncgene = [x for x in keys_tgenes if self.nc_genes[index] in x]
                ncgenes_all = ncgenes_all + ncgene
            ncgenes_values = np.array(list(dict.fromkeys(ncgenes_all)))                     # remove duplicate entries
        else:
            raise TypeError("Negative control genes must be supplied as str or list of strings.")
        return ncgenes_values


    def prepare_nc_model(self, keys_tgenes, strong, sens_lethality, sens_recovery, nc_values, cell):
        '''
        Function to define a negative control model
        :param keys_tgenes: Numpy array with names of all triple gene combinations
        :param strong: Numpy array containing strong scores
        :param sens_lethality: Numpy array containing sensitive lethality scores
        :param sens_recovery: Numpy array containing sensitive recovery scores
        :param nc_values: np array containing names of triple gene combinations which include at least one nc gene
        :param cell: Integer indicating the cell line
        :return: 3 np arrays containing the scores of the negative control combinations
                 (for strong, sensitive lethality and sensitive recovery scores)
        '''
        strong_dic = dict(zip(keys_tgenes, strong[:, cell]))         # create new dictionary for score values
        strong_nc = np.array([strong_dic[k] for k in nc_values])     # obtain only values whose keys are nc combinations

        if len(self.pc_genes) > 0:
            lethality_dic = dict(zip(keys_tgenes, sens_lethality[:, cell]))
            recovery_dic = dict(zip(keys_tgenes, sens_recovery[:, cell]))
            lethality_nc = np.array([lethality_dic[k] for k in nc_values])
            lethality_nc = lethality_nc[~np.isnan(lethality_nc)]                                # drop nan
            recovery_nc = np.array([recovery_dic[k] for k in nc_values])
            recovery_nc = recovery_nc[~np.isnan(recovery_nc)]                                   # drop nan
        else:
            lethality_nc = np.zeros(np.shape(strong_nc))
            recovery_nc = np.zeros(np.shape(strong_nc))

        return strong_nc, lethality_nc, recovery_nc


    def gmm(self, nc_input, title, cell):
        '''
        Function to find a Gaussian mixture model (GMM) which best describes the distribution of negative control scores.
        sklearn.mixture and GMMs with 1-5 Gaussians are used.
        :param nc_input: Numpy array containing the scores of the negative control combinations
        :param title: 'Strong', 'lethality' or 'recovery'
        :param cell: Name of cell line as string
        :return: GMM which was considered as best solution by the Bayesian information criterion.
        '''
        nc = np.reshape(nc_input, (-1, 1))                                      # reshaping required for fit

        # fit Gaussian Mixture Models with 1 to max. 5 components, respectively, and find the best fit
        if np.size(nc) < 5:
            number = np.arange(1, np.size(nc) + 1)                              # array with length of nc array
        else:
            number = np.arange(1, 6)                                            # array with 5 elements
        models = [None for i in range(len(number))]                             # list with 5 entries "None"

        for i in range(len(number)):
            models[i] = GaussianMixture(n_components=number[i], tol=0.0001, max_iter=100, n_init=3,
                                        init_params='kmeans', verbose=0).fit(nc)

        bic = [m.bic(nc) for m in models]                                       # Bayesian information criterion
        best_model = models[bic.index(min(bic))]                                # Minimum of BIC defines best model

        # Print information about current status in terminal
        if models[bic.index(min(bic))].converged_ is True and self.verbose:
            print("Gaussian Mixture Model converged.")
        else:
            print("No convergence was reached for Gaussian Mixture Model.")
        txt = "According to the Bayesian information criterion, {0} steps were required" \
              " by the best fit of EM to converge."
        if self.verbose:
            print(txt.format(models[bic.index(min(bic))].n_iter_))
        if np.argmin(bic) == 0:
            print("Best GMM consists of 1 Gaussian.")
        else:
            print("Best GMM consists of", np.argmin(bic)+1, "Gaussians.")

        if self.plotting:
            # Plot figures showing the result of fitting GMM to the data
            fig = plt.figure(figsize=(12, 8))
            fig.subplots_adjust(left=0.12, right=0.9, bottom=0.21, top=0.9, wspace=0.5)

            # plot 1: data + best-fit mixture
            ax = fig.add_subplot(121)                                               # subplot 1 (left)
            m_best = models[np.argmin(bic)]
            lin = np.linspace(min(nc) - 0.1, max(nc) + 0.1, num=100)
            logprob = m_best.score_samples(lin.reshape(-1, 1))
            responsibilities = m_best.predict_proba(lin.reshape(-1, 1))
            pdf = np.exp(logprob)                                                   # probability density
            pdf_individual = responsibilities * pdf[:, np.newaxis]
            ax.hist(nc, 30, density=True, histtype='stepfilled', color='#009FE3', alpha=0.4)
            if np.argmin(bic) == 0:                                                 # GMM with only one component
                ax.plot(lin, pdf, '--k')                                            # display model without legend
            else:                                                                   # GMMs with 2 - 5 components
                def legend(ax):                                                     # legend without redundant labels
                    handles, label = ax.get_legend_handles_labels()
                    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, label)) if l not in label[:i]]
                    ax.legend(*zip(*unique), loc=2, fontsize=14)
                ax.plot(lin, pdf, '--k', label='Gaussian mixture')
                ax.plot(lin, pdf_individual, ':k', label='individual Gaussians')
                legend(ax)
            ax.text(0.0, 1.08, "Best-fit GMM " + title + cell, ha='left', va='top', transform=ax.transAxes, fontsize=20)
            ax.set_xlabel('score nc gene combinations', fontsize=18)
            ax.set_ylabel('count', fontsize=20)
            ax.tick_params(axis='both', labelsize=14)

            # plot 2: BIC
            ax = fig.add_subplot(122)                                               # subplot 2 (right)
            ax.plot(number, bic, color='#009FE3', label='BIC')
            ax.set_xticks(np.arange(min(number), max(number) + 1, 1.0))             # define step size of x ticks
            ax.text(0.0, 1.08, 'Number of Gaussians ' + cell, ha='left', va='top', transform=ax.transAxes, fontsize=20)
            ax.set_xlabel('components', fontsize=20)
            ax.set_ylabel('information criterion', fontsize=20)
            ax.tick_params(axis='both', labelsize=14)
            ax.legend(loc=2, fontsize=16)                                           # legend on the top left

            plt.savefig(os.path.join(self.output_dir, "Output", "GMM", ("GMM_" + title + "_" + cell + ".png")))
            plt.close()
        return best_model


    def p_value(self, nc_model, score):
        '''
        Function for calculating P-values for each individual score. A P-value describes the probability of the test
        result given an actual distribution of negative controls (approximated by the Gaussian mixture model).
        p_values = 1 - (weight * (norm(mean, stDev).cdf(score)))
        :param nc_model: Gaussian mixture models which best describes the distribution of nc genes
        :param score: array of all 'strong', 'lethality' or 'recovery' scores, respectively
        :return: Numpy array of P-values in same order as scores
        '''
        n = np.size(nc_model.means_)                                            # get total number of elements in array
        added_elements = 0
        for i in range(0, n, 1):
            element = nc_model.weights_[i] * (norm(nc_model.means_[i], np.sqrt(nc_model.covariances_[i])).sf(score))
            added_elements = added_elements + element
        pvalue = added_elements
        return np.ravel(pvalue)


    def fdr(self, pvalue):
        '''
        Function for calculating adjusted P-values (FDR) by using the method of Benjamini and Hochberg.
        P-values are sorted before calculating the FDRs.
        :param pvalue: Numpy array of P-values of one score
        :return: Numpy array of adjusted P-values of one score
        '''
        pvalue_sort = np.sort(pvalue, axis=None)                                # returns the flattened array, sorted
        pvalue_argsort = np.argsort(pvalue, axis=None)                          # returns original indices of array
        mask = np.isfinite(pvalue_sort)                                         # mask for pvalues which are not NAN

        # Benjamini/Hochberg (non-negative)
        fdr_bh = multi.multipletests(pvalue_sort[mask], alpha=0.05, method='fdr_bh', is_sorted=True)

        fdr_bh_sorted = np.empty_like(pvalue)  # empty array to put back in the NAN
        fdr_bh_sorted_out = np.empty_like(pvalue)  # empty array to return
        fdr_bh_sorted.fill(np.nan)
        fdr_bh_sorted[mask] = fdr_bh[1]

        fdr_bh_sorted_out[pvalue_argsort] = fdr_bh_sorted  # return fdr to original order
        return fdr_bh_sorted_out


    def dataframe(self, score, pvalue, fdr, cell_list, line, names):
        '''
        Function for creating an empty output data frame.
        :param score: array of all 'strong', 'lethality' or 'recovery' scores, respectively
        :param pvalue: Numpy array of P-values
        :param fdr: Numpy array of adjusted P-values
        :param cell_list: String or list of cell lines included in CRISPR screen
        :param line: Integer indicating the cell line
        :param names: Data frame where triple genes are mapped to their single genes
        :return: Data frame with three columns (Score, P-value, Adj_P-value) for the first cell line
        '''
        df = pd.DataFrame(data=np.column_stack((score, pvalue, fdr)), index=names.index.tolist(),
                          columns=['Score_' + cell_list[line],
                                   'P-value_' + cell_list[line],
                                   'Adj_P-value_' + cell_list[line]])
        return df


    def add_columns(self, df, score, pvalue, fdr, cell_list, line):
        '''
        Function for adding three columns to the output data frame.
        :param df: Data frame with three columns (Score, P-value, Adj_P-value) for every cell line
        :param score: array of all 'strong', 'lethality' or 'recovery' scores, respectively
        :param pvalue: Numpy array of P-values
        :param fdr: Numpy array of adjusted P-values
        :param cell_list: String or list of cell lines included in CRISPR screen
        :param line: Integer indicating the cell line
        :return: Data frame with three new columns for every additional cell line
        '''
        df['Score_' + cell_list[line]] = score[:, line].tolist()
        df['P-value_' + cell_list[line]] = pvalue.tolist()
        df['Adj_P-value_' + cell_list[line]] = fdr.tolist()
        return df


    def get_top10_df(self, df, cellpos, cellline):
        '''
        Function for creating a data frame for the 10 gene triples with the strongest effects.
        :param df: Output data frame of the strong score
        :param cellpos: Integer indicating the cell line
        :param cellline: Name of the cell line as string
        :return: Data frame with 10 rows for the most interesting triple combinations with a lethal, recovering,
                 or no effect (directly saved in output folder)
        '''
        if df.iloc[:, cellpos * 4].str.contains('lethality').any():
            lethality_df = df[df.iloc[:, cellpos * 4] == 'lethality'].iloc[:, cellpos * 4:cellpos * 4 + 4]
            top10_lethal = lethality_df.sort_values(by=[lethality_df.columns[1]], ascending=False).iloc[0:10, :]
            top10_lethal.to_csv(os.path.join(self.output_dir, "Output", ("top10_lethal_" + cellline + ".csv")))
            top10_lethal.name = 'top10_lethal_' + cellline
            self.scatter(top10_lethal)

        if df.iloc[:, cellpos * 4].str.contains('recovery').any():
            recovery_df = df[df.iloc[:, cellpos * 4] == 'recovery'].iloc[:, cellpos * 4:cellpos * 4 + 4]
            top10_recovery = recovery_df.sort_values(by=[recovery_df.columns[1]], ascending=False).iloc[0:10, :]
            top10_recovery.to_csv(os.path.join(self.output_dir, "Output", ("top10_recovery_" + cellline + ".csv")))
            top10_recovery.name = 'top10_recovery_' + cellline
            self.scatter(top10_recovery)

        if df.iloc[:, cellpos * 4].str.contains('no').any():
            noeffect_df = df[df.iloc[:, cellpos * 4] == 'no effect'].iloc[:, cellpos * 4:cellpos * 4 + 4]
            top10_noeffect = noeffect_df.sort_values(by=[noeffect_df.columns[1]], ascending=False).iloc[0:10, :]
            top10_noeffect.to_csv(os.path.join(self.output_dir, "Output", ("top10_noeffect_" + cellline + ".csv")))
            top10_noeffect.name = 'top10_noeffect_' + cellline
            self.scatter(top10_noeffect)


    def scoring(self):
        '''
        Main function of the class Score.
        All the other functions described above are called from this function.
        :return: - 3 final data frames, one for each score (columns for effect, score, P-value, and FDR per cell line)
                 - Cell lines included in CRISPR screen (as integer or list)
        '''
        # Create new directory for outputs
        if not os.path.exists(os.path.join(self.output_dir, "Output")):
            os.makedirs(os.path.join(self.output_dir, "Output"))

        # Create new directory for results and plots of Gaussian Mixture Model
        if not os.path.exists(os.path.join(self.output_dir, "Output", "GMM")):
            os.makedirs(os.path.join(self.output_dir, "Output", "GMM"))

        # Create new directory for plots
        if not os.path.exists(os.path.join(self.output_dir, "Output", "Plots")):
            os.makedirs(os.path.join(self.output_dir, "Output", "Plots"))

        # Modify structure of output data from of Model1 (create informative dictionaries)
        mod_y = self.new_dic(self.model.y, self.model.hash_y)
        mod_s = self.new_dic(self.model.s, self.model.hash_s)

        # Create data frame where triple genes are mapped to their single genes
        Tgenes2Sgenes = self.df_triple(mod_s)
        keys_tgenes = Tgenes2Sgenes.index.to_numpy()                                      # convert rownames to np array

        # Get information about which cell lines are present in data
        cell_list = self.cell_info(mod_y)

        # Create arrays with y values in correct order
        yg = self.get_y(Tgenes2Sgenes.iloc[:, 0].values, mod_y, cell_list)                # use columns of Tgenes2Sgenes
        yh = self.get_y(Tgenes2Sgenes.iloc[:, 1].values, mod_y, cell_list)                # and values from mod_y array
        yf = self.get_y(Tgenes2Sgenes.iloc[:, 2].values, mod_y, cell_list)
        sghf = self.get_s(Tgenes2Sgenes.index,mod_s,cell_list)

        # Compute scores
        strong, sensitive_lethality, sensitive_recovery = self.calculate_scores(yg, yh, yf, sghf, mod_y)

        # Create np array with values of negative control genes
        ncgenes_values = self.nc_values(keys_tgenes)

        # Create nc_model per cell line to calculate p-value and FDR per cell line
        for cell_line in range(len(cell_list)):

            # Obtain negative control gene values within each score (as input for Gaussian Mixture Model)
            strong_nc, lethality_nc, recovery_nc = self.prepare_nc_model(keys_tgenes, strong, sensitive_lethality,
                                                                         sensitive_recovery, ncgenes_values, cell_line)

            # Check distribution of negative control gene score values before running EM algorithm
            if self.plotting:
                self.nc_distribution(strong_nc, lethality_nc, recovery_nc, cell_list[cell_line])

            # Running expectation maximization algorithm from Gaussian Mixture on the three scores
            strong_nc_model = self.gmm(strong_nc, "strong ", cell_list[cell_line])
            if len(self.pc_genes) > 0:
                sensitive_lethality_nc_model = self.gmm(lethality_nc, "lethality ", cell_list[cell_line])
                sensitive_recovery_nc_model = self.gmm(recovery_nc, "recovery ", cell_list[cell_line])

            # Calculate P-values (null hypothesis: strong_nc_model > strong)
            pvalue_strong = self.p_value(strong_nc_model, strong[:, cell_line])
            if len(self.pc_genes) > 0:
                pvalue_lethality = self.p_value(sensitive_lethality_nc_model, sensitive_lethality[:, cell_line])
                pvalue_recovery = self.p_value(sensitive_recovery_nc_model, sensitive_recovery[:, cell_line])

            # Calculate Adjusted P-values (FDR)
            fdr_strong = self.fdr(pvalue_strong)
            if len(self.pc_genes) > 0:
                fdr_lethality = self.fdr(pvalue_lethality)
                fdr_recovery = self.fdr(pvalue_recovery)

            # Extract information whether strong combinatorial effect is lethal or recovering
            effect_type = []
            for i in range(len(strong[:, cell_line])):                                  # iterate through rows of strong
                if strong[i, cell_line] > 0 and sghf[i, cell_line] >= 0:
                    effect_type.append('recovery')
                elif strong[i, cell_line] > 0 and sghf[i, cell_line] < 0:
                    effect_type.append('lethality')
                else:
                    effect_type.append('no effect')
            strong_type = np.array(effect_type)

            # Create output df (with initially three columns, and 3 extra columns per additional cell line)
            if cell_line == 0:
                strong_df = self.dataframe(strong[:, cell_line], pvalue_strong, fdr_strong,
                                           cell_list, cell_line, Tgenes2Sgenes)
                strong_df.insert(loc=0, column='Interaction_' + cell_list[cell_line], value=strong_type)
                if len(self.pc_genes) > 0:
                    sensitive_lethality_df = self.dataframe(sensitive_lethality[:, cell_line], pvalue_lethality,
                                                            fdr_lethality, cell_list, cell_line, Tgenes2Sgenes)
                    sensitive_recovery_df = self.dataframe(sensitive_recovery[:, cell_line], pvalue_recovery,
                                                           fdr_recovery, cell_list, cell_line, Tgenes2Sgenes)
                else:  # create empty data frames (required for order of elements returned by function "scoring()"
                    data = {'Score_' + cell_list[cell_line]: [0],
                            'P-value_' + cell_list[cell_line]: [0],
                            'Adj_P-value_' + cell_list[cell_line]: [0]}
                    sensitive_lethality_df = pd.DataFrame(data, index=['no data available'])
                    sensitive_recovery_df = pd.DataFrame(data, index=['no data available'])

            else:
                strong_df.insert(loc=cell_line*3+1, column='Interaction_' + cell_list[cell_line], value=strong_type)
                strong_df = self.add_columns(strong_df, strong, pvalue_strong, fdr_strong, cell_list, cell_line)
                if len(self.pc_genes) > 0:
                    sensitive_lethality_df = self.add_columns(sensitive_lethality_df, sensitive_lethality,
                                                              pvalue_lethality, fdr_lethality, cell_list, cell_line)
                    sensitive_recovery_df = self.add_columns(sensitive_recovery_df, sensitive_recovery,
                                                             pvalue_recovery, fdr_recovery, cell_list, cell_line)

            # Extract important information from output dfs
            self.get_top10_df(strong_df, cell_line, cell_list[cell_line])

        # Save output data frames
        strong_df.to_csv(os.path.join(self.output_dir, "Output", "strong_df.csv"))
        if len(self.pc_genes) > 0:
            sensitive_lethality_df.to_csv(os.path.join(self.output_dir, "Output", "sensitive_lethality_df.csv"))
            sensitive_recovery_df.to_csv(os.path.join(self.output_dir, "Output", "sensitive_recovery_df.csv"))

        print("Scores were successfully calculated.")
        return strong_df, sensitive_lethality_df, sensitive_recovery_df, cell_list


    # The second part of the class Score contains functions for creating visualizations
    def convergence(self, mae):
        '''
        Function to create a line plot showing the convergence of the CAVI algorithm.
        If model is divergent, alternative parameters for the priors should be selected until convergence is achieved.
        :param mae: Mean absolute error after every iteration (~ convergence rate).
        :return: Figure (saved in output folder)
        '''
        plt.figure(figsize=(12, 8))
        plt.plot(mae, color='#009FE3')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.title('Convergence of interference', fontsize=22)
        plt.xlabel('Iteration steps', fontsize=20)
        plt.ylabel('Mean absolute error', fontsize=20)
        plt.savefig(os.path.join(self.output_dir, "Output", "GMM", "Convergence.png"))
        plt.close()


    def nc_distribution(self, strong_nc, sens_lethality_nc, sens_recovery_nc, cell):
        '''
        Function to depict the distribution of negative control genes (should be around 0 = no effect).
        :param strong_nc: Numpy array of strong scores containing negative control genes.
        :param sens_lethality_nc: Numpy array of sensitive lethality scores containing negative control genes.
        :param sens_recovery_nc: Numpy array of sensitive recovery scores containing negative control genes.
        :param cell: Name of cell line as string.
        :return: Figure showing three box plots (saved in output folder)
        '''
        plt.figure(figsize=(12, 8))

        if len(self.pc_genes) > 0:
            box = plt.boxplot([strong_nc, sens_lethality_nc, sens_recovery_nc], patch_artist=True, showmeans=True)
            plt.setp(box['means'], marker='D')
            c = ['#0047B9', '#FF5C00', '#00CC00']
            c2 = ['#0047B9', '#0047B9', '#FF5C00', '#FF5C00', '#00CC00', '#00CC00']
            c3 = ['#009FE3', '#C61826', 'darkgreen']
            ticks = [1, 2, 3]
            labels = ['strong', 'sensitive\nlethality', 'sensitive\nrecovery']
        else:
            box = plt.boxplot([strong_nc], patch_artist=True, showmeans=True)
            plt.setp(box['means'], marker='D')
            c = ['#0047B9']
            c2 = ['#0047B9', '#0047B9']
            c3 = ['#009FE3']
            ticks = [1]
            labels = ['strong']

        # Define colors of box plots
        for patch, color in zip(box['boxes'], c):
            patch.set_color(color)
            patch.set_facecolor(color)
        for patch, color in zip(box['fliers'], c):
            patch.set_color(color)
            patch.set_markeredgecolor(color)
        for patch, color in zip(box['whiskers'], c2):
            patch.set_color(color)
        for patch, color in zip(box['caps'], c2):
            patch.set_color(color)
        for patch, color in zip(box['means'], c3):
            patch.set_markeredgecolor(color)
            patch.set_markerfacecolor(color)
        for patch, color in zip(box['medians'], c3):
            patch.set_color(color)

        plt.xticks(ticks, labels, fontsize=20)
        plt.yticks(fontsize=14)
        plt.ylabel('scores', fontsize=20)
        title = 'Distribution of negative control genes ' + cell
        plt.title(title, fontsize=20)
        plt.legend([box['medians'][0], box['means'][0]], ['median', 'mean'], fontsize=16)
        title = title.replace(" ", "_")
        plt.savefig(os.path.join(self.output_dir, "Output", "GMM", (title+".png")))
        plt.close()


    def LFC_distribution(self):
        '''
        Function to show the distribution of log2 fold changes.
        :return: Box plot of LFCs (saved in output folder)
        '''
        plt.figure(figsize=(12, 8))
        c = "#0047B9"
        box = plt.boxplot(self.model.D, patch_artist=True, showmeans=True,
                          boxprops=dict(facecolor=c, color=c),
                          capprops=dict(color=c),
                          whiskerprops=dict(color=c),
                          flierprops=dict(color=c, markeredgecolor=c),
                          medianprops=dict(color='#009FE3'),
                          meanprops=dict(markeredgecolor='#009FE3', markerfacecolor='#009FE3', marker='D'))
        ticks = list(range(len(self.score_result[3]) + 1))                             # create 1 tick per cell line
        ticks.pop(0)                                                                   # delete first element which is 0
        plt.xticks(ticks, self.score_result[3], fontsize=18)
        plt.ylabel('$Log_2$ fold changes', fontsize=20)
        plt.yticks(fontsize=14)
        title = 'Distribution of LFC'
        plt.title(title, fontsize=20)
        plt.legend([box['medians'][0], box['means'][0]], ['median', 'mean'], fontsize=16)
        plt.savefig(os.path.join(self.output_dir, "Output", "Plots", (title + ".png")))
        plt.close()


    def scatter(self, top10df):
        '''
        Function to create a scatter plot showing the 10 gene combinations with the strongest effect.
        :param top10df: Data frame with 10 rows for the most lethal or recovering triple combinations
        :return: Scatter plot (saved in output folder)
        '''
        effect = top10df.name.split('_')[1]                           # effect = 'lethality', 'recovery', or 'no effect'
        cell = top10df.name.split('_')[2]                                               # get cell line
        y = top10df.iloc[:, 1]
        x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        x = x[0:np.size(y)]                                                             # adjust length of x by y
        if effect == 'lethal':                                                          # orange color palette
            colors = ['#A83C00', '#BF4500', '#F25600', '#FF5C00', '#FF9053',
                      '#FF9D66', '#FFA675', '#FFBE99', '#FFB78F', '#FFDECC']
        elif effect == 'recovery':                                                      # green color palette
            colors = ['#008600', '#009E00', '#00B000', '#00C000', '#00C800',
                      '#05FF05', '#3FFF3F', '#5DFF5D', '#89FF89', '#A7FFA7']
        else:                                                                           # blue color palette (no effect)
            colors = ['#006D9E', '#0077AA', '#0094D6', '#009FE3', '#00AAF6',
                      '#47C6FF', '#55CCFF', '#75D4FF', '#8EDDFF', '#A7E4FF']
        plt.figure(figsize=(12, 8))
        plt.scatter(x, y, s=350, c=colors[0:np.size(y)])
        plt.xticks(x, top10df.index, rotation=45, fontsize=14)                          # set text labels and properties
        plt.yticks(fontsize=14)
        plt.margins(0.1)
        plt.subplots_adjust(bottom=0.35)
        plt.ylabel('score', fontsize=20)
        if effect == 'noeffect':
            title = 'Top ' + str(np.size(y)) + ' triple gene interactions with no effect ' + cell
        else:
            title = 'Top ' + str(np.size(y)) + ' ' + effect + ' triple gene interactions ' + cell
        plt.title(title, fontsize=20)
        title = title.replace(" ", "_")
        plt.savefig(os.path.join(self.output_dir, "Output", "Plots", (title+".png")))
        plt.close()


    def volcano_score(self, df, effect, col=None):
        '''
        Function to create a volcano plot showing the distribution and highlighting interesting calculated scores.
        :param df: Final data frame for score with the columns 'effect', 'score', 'P-value', and 'FDR' per cell line
        :param effect: String describing effect ('strong', 'sensitive lethality', or 'sensitive recovery')
        :return: Volcano plot showing distribution of scores (saved in output folder)
        '''
        col = [] if col is None else col
        if effect == 'strong':
            x = df.iloc[:, 1]
            y = (df.iloc[:, 1]) ** 2
            for i in (range(len(x))):
                col.append(np.where(df.iloc[i, 0] == 'lethality', '#FF5C00',
                                    np.where(df.iloc[i, 0] == 'recovery', '#00CC00', '#009FE3')))
        elif effect == 'sensitive lethality':
            x = df.iloc[:, 0]
            y = (df.iloc[:, 0]) ** 2
            for i in (range(len(x))):
                col.append(np.where(pd.isna(df.iloc[i, 0]), '#009FE3', '#FF5C00'))
        else:
            x = df.iloc[:, 0]
            y = (df.iloc[:, 0]) ** 2
            for i in (range(len(x))):
                col.append(np.where(pd.isna(df.iloc[i, 0]), '#009FE3', '#00CC00'))
        color = np.array(col)
        plt.figure(figsize=(12, 8))
        title = 'Distribution of ' + effect + ' scores'
        plt.title(title, fontsize=21)
        plt.xlabel('score', fontsize=20)
        plt.xticks(fontsize=14)
        plt.ylabel('$score^2$', fontsize=20)
        plt.yticks(fontsize=14)
        plt.scatter(x, y, s=120, c=color, alpha=0.6, edgecolors='none')                     # alpha: transparency
        a = mpatches.Patch(color='#FF5C00', label='lethality')
        b = mpatches.Patch(color='#00CC00', label='recovery')
        c = mpatches.Patch(color='#009FE3', label='no effect')
        plt.legend(handles=[a, b, c], fontsize=16)
        title = title.replace(" ", "_")
        plt.savefig(os.path.join(self.output_dir, "Output", "Plots", (title+".png")))
        plt.close()

    def volcano_LFC(self, df, col=None):
        '''
        Function to create a volcano plot showing the mean combinatorial gene effect vs -log10(P-value).
        In this way the effect sizes of the triple gene combinations can be compared.
        :param df: Final data frame for score with the columns 'effect', 'score', 'P-value', and 'FDR' per cell line
        :return: Volcano plot showing LFCs (saved in output folder)
        '''
        col = [] if col is None else col
        x = self.model.s.mean(axis=1)                                                               # mean log2 FC
        y = np.abs(np.log10(df.iloc[:, 2]))                                                         # -log10(P-value)
        for i in (range(len(x))):
            col.append(np.where(df.iloc[i, 0] == 'lethality', '#FF5C00',
                                np.where(df.iloc[i, 0] == 'recovery', '#00CC00', '#009FE3')))
        color = np.array(col)
        plt.figure(figsize=(12, 8))
        plt.scatter(x, y, s=120, c=color, alpha=0.6, edgecolors='none')                     # alpha: transparency
        title = 'Volcano strong combinatorial effect'
        plt.title(title, fontsize=21)
        plt.xlabel('Mean combinatorial effect', fontsize=20)
        plt.xticks(fontsize=14)
        plt.ylabel('$-log_{10}$(P-value)', fontsize=20)
        plt.yticks(fontsize=14)
        a = mpatches.Patch(color='#FF5C00', label='lethality')
        b = mpatches.Patch(color='#00CC00', label='recovery')
        c = mpatches.Patch(color='#009FE3', label='no effect')
        plt.legend(handles=[a, b, c], fontsize=14)
        plt.hlines(1.3, min(x), max(x), color='grey', linestyle='dashed')                   # plotting line at p=0.05
        plt.text(max(x), 1.3, 'p<0.05', fontsize=11, color='grey',
                 va='center_baseline', ha='center', backgroundcolor='w')

        # Create labels for the 2 most important lethal as well as recovering interactions
        pvalue = df.iloc[:, 2].to_numpy()                                            # extract P-values in correct order
        top2_lethality = df.sort_values(by=[df.columns[0], df.columns[1]], ascending=[True, False]).iloc[0:2, :]
        top2_recovery = df.sort_values(by=[df.columns[0], df.columns[1]], ascending=[False, False]).iloc[0:2, :]
        top2_lethality = top2_lethality[top2_lethality.iloc[:, 0] == 'lethality']    # keep only lethal interactions
        top2_recovery = top2_recovery[top2_recovery.iloc[:, 0] == 'recovery']        # keep only recovery interactions
        top2_lethality_np = top2_lethality.iloc[:, 2].to_numpy()                     # convert P-values to np array
        top2_recovery_np = top2_recovery.iloc[:, 2].to_numpy()

        pos_list_lethal = []
        pos_list_recovery = []
        for pos, value in enumerate(pvalue):
            for val_l in top2_lethality_np:
                if value == val_l:
                    pos_list_lethal.append(pos)
            for val_r in top2_recovery_np:
                if value == val_r:
                    pos_list_recovery.append(pos)

        # extract rownames for annotation of data points
        rownames = list(df.index)
        if pos_list_lethal == pos_list_recovery:                                    # in case there is only 'no effect'
            for pos in pos_list_lethal:
                plt.annotate(rownames[pos], xy=(x[pos], y[pos]), fontsize=12,
                             arrowprops=dict(fc='#009FE3', width=0.9, shrink=0.05),
                             horizontalalignment='center', verticalalignment='top',
                             bbox=dict(fc='#FF5C00', alpha=0.05, ec='#009FE3'))
        else:
            for pos in pos_list_lethal:
                plt.annotate(rownames[pos], xy=(x[pos], y[pos]), fontsize=12,
                             arrowprops=dict(fc='#FF5C00', width=0.9, shrink=0.05),
                             horizontalalignment='left', verticalalignment='top',
                             bbox=dict(fc='#FF5C00', alpha=0.05, ec='#FF5C00'))

            for pos in pos_list_recovery:
                plt.annotate(rownames[pos], xy=(x[pos], y[pos]), fontsize=12,
                             arrowprops=dict(fc='#00CC00', width=0.9, shrink=0.05),
                             horizontalalignment='right', verticalalignment='top',
                             bbox=dict(fc='#00CC00', alpha=0.05, ec='#00CC00'))

        title = title.replace(" ", "_")
        plt.savefig(os.path.join(self.output_dir, "Output", "Plots", (title+".png")))
        plt.close()


    def plotting_function(self):
        '''
        Main function to call all functions that create a figure.
        '''
        # Mean absolute error
        self.convergence(self.model.MAE)

        # Distribution log2 fold changes of triple gene interactions
        self.LFC_distribution()

        # Volcano plot LFC
        self.volcano_LFC(self.score_result[0])

        # Volcano plot with scores
        self.volcano_score(self.score_result[0], 'strong')
        if len(self.pc_genes) > 0:
            self.volcano_score(self.score_result[1], 'sensitive lethality')
            self.volcano_score(self.score_result[2], 'sensitive recovery')

        print("Plots were successfully drawn.\n"
              "Figures and resulting data frames can be seen in folder 'Output'.")
