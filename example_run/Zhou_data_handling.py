### Explanation how the raw data from Zhou were converted into the inputs for Gemi3 (LFC_table and annotation_table)

'''
1. Download all twelve SRA-files from the GEO accession viewer (GSM4664050 - GSM4664061)
   (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154112)
   Example for GSM4664050: https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR12185790

2. Convert all SRA files into FASTQ files by using the fastq-dump command from the SRA Toolkit on each SRA file
   (for example: sratoolkit.2.11.1-win64)

3. Merge all FASTQ files into one large FASTQ file (~33 GB)

4. Install CombiPIPE (https://github.com/AWHKU/CombiPIPE)

5. Run CombiPIPE with the FASTQ file, sample_info.csv and barcode_list.csv as input

6. Now, we need the eight csv files in the folder outputs_newBCanalyzer
   (BCcount_SA160.fastq_TCCAGGCAT.csv - BCcount_SA166.fastq_GATGCGATTGCG.csv)
   The column 'BCcounts' contains the raw counts of every gene triple.
   As the keys and BCcounts are in different orders in every file, they need to be sorted and merged into one file.
   Therefore, the columns 'key' and 'BCcounts' need to be renamed in each file so that they are not redundant.
   We create a new file 'counts_empty.csv' that has one column 'key' with strings from 1_1_1 to 32_32_32 (32768 rows).
'''

import pandas as pd
import itertools as itt

# Import data frames

counts_empty = pd.DataFrame({"key":[f"{i[0]}_{i[1]}_{i[2]}" for i in itt.product(range(1,33),repeat=3)]})
counts1 = pd.read_csv("BCcount_SA160.fastq_TCCAGGCAT.csv", sep=",",skiprows=7)
counts1.columns=['BC3', 'BC2', 'BC1', 'key_160_1', 'sgRNAs', 'genes', 'BCcounts_160_1', 'log2CPM']
counts1 = counts1.loc[:,['key_160_1', 'BCcounts_160_1']]

counts2 = pd.read_csv("BCcount_SA160.fastq_TCCAGGCATGGA.csv", sep=",",skiprows=7)
counts2.columns=['BC3', 'BC2', 'BC1', 'key_160_2', 'sgRNAs', 'genes', 'BCcounts_160_2', 'log2CPM']
counts2 = counts2.loc[:,['key_160_2', 'BCcounts_160_2']]

counts3 = pd.read_csv("BCcount_SA162.fastq_GATCAATGT.csv", sep=",",skiprows=7)
counts3.columns=['BC3', 'BC2', 'BC1', 'key_162_1', 'sgRNAs', 'genes', 'BCcounts_162_1', 'log2CPM']
counts3 = counts3.loc[:,['key_162_1', 'BCcounts_162_1']]

counts4 = pd.read_csv("BCcount_SA162.fastq_GATCAATGTTCG.csv", sep=",",skiprows=7)
counts4.columns=['BC3', 'BC2', 'BC1', 'key_162_2', 'sgRNAs', 'genes', 'BCcounts_162_2', 'log2CPM']
counts4 = counts4.loc[:,['key_162_2', 'BCcounts_162_2']]


counts5 = pd.read_csv("BCcount_SA164.fastq_TCGTTCCTG.csv", sep=",",skiprows=7)
counts5.columns=['BC3', 'BC2', 'BC1', 'key_164_1', 'sgRNAs', 'genes', 'BCcounts_164_1', 'log2CPM']
counts5 = counts5.loc[:,['key_164_1', 'BCcounts_164_1']]


counts6 = pd.read_csv("BCcount_SA164.fastq_TCGTTCCTGGGA.csv", sep=",",skiprows=7)
counts6.columns=['BC3', 'BC2', 'BC1', 'key_164_2', 'sgRNAs', 'genes', 'BCcounts_164_2', 'log2CPM']
counts6 = counts6.loc[:,['key_164_2', 'BCcounts_164_2']]


counts7 = pd.read_csv("BCcount_SA166.fastq_GATGCGATT.csv", sep=",",skiprows=7)
counts7.columns=['BC3', 'BC2', 'BC1', 'key_166_1', 'sgRNAs', 'genes', 'BCcounts_166_1', 'log2CPM']
counts7 = counts7.loc[:,['key_166_1', 'BCcounts_166_1']]


counts8 = pd.read_csv("BCcount_SA166.fastq_GATGCGATTGCG.csv", sep=",",skiprows=7)
counts8.columns=['BC3', 'BC2', 'BC1', 'key_166_2', 'sgRNAs', 'genes', 'BCcounts_166_2', 'log2CPM']
counts8 = counts8.loc[:,['key_166_2', 'BCcounts_166_2']]

# Merge data frames by restoring the correct order of keys
counts_df = counts_empty.merge(counts1, left_on='key', right_on='key_160_1', how='left')
counts_df = counts_df.merge(counts2, left_on='key', right_on='key_160_2', how='left')
counts_df = counts_df.merge(counts3, left_on='key', right_on='key_162_1', how='left')
counts_df = counts_df.merge(counts4, left_on='key', right_on='key_162_2', how='left')
counts_df = counts_df.merge(counts5, left_on='key', right_on='key_164_1', how='left')
counts_df = counts_df.merge(counts6, left_on='key', right_on='key_164_2', how='left')
counts_df = counts_df.merge(counts7, left_on='key', right_on='key_166_1', how='left')
counts_df = counts_df.merge(counts8, left_on='key', right_on='key_166_2', how='left')

# We keep only the BCcounts columns because we only need one key column
counts_ordered = pd.DataFrame(counts_df, columns=['key', 'BCcounts_160_1', 'BCcounts_160_2', 'BCcounts_162_1', 'BCcounts_162_2',
                                                  'BCcounts_164_1', 'BCcounts_164_2', 'BCcounts_166_1', 'BCcounts_166_2'])

# The columns need to be renamed in order to run 'calc_observed_LFC' in g3datahandling.py
counts_renamed = counts_ordered.reindex(columns=['key', 'Cell_D15.rep1', 'Cell_D15.rep2', 'Cell_D26.rep1', 'Cell_D26.rep2',
                                                 'Cell_D15.rep3', 'Cell_D15.rep4', 'Cell_D26.rep3', 'Cell_D26.rep4'])
counts_renamed.to_csv("/home/l457h/counts_ordered.csv", sep="\t", index=False)

'''
7. In order to run run 'calc_observed_LFC' in g3datahandling.py, three inputs are required:
   - counts_ordered.csv as input 'raw_data'
   - ETP = ["Cell_D15.rep1","Cell_D15.rep2","Cell_D15.rep3","Cell_D15.rep4"]
   - replicates_2D_zhou.txt as input 'replicates':
     colname	    samplename	replicate
     Cell_D15.rep1	Cell_D15	rep1
     Cell_D15.rep2	Cell_D15	rep2
     Cell_D26.rep1	Cell_D26	rep1
     Cell_D26.rep2	Cell_D26	rep2
     Cell_D15.rep3	Cell_D15	rep3
     Cell_D15.rep4	Cell_D15	rep4
     Cell_D26.rep3	Cell_D26	rep3
     Cell_D26.rep4	Cell_D26	rep4

8. The output of 'calc_observed_LFC' can be saved as LFC_raw_counts.csv

9. Create the three input files for Gemi3:
   - LFC_table.csv (table containing keys and log fold changes, all empty/nan rows need to be deleted)
   - samples = ["OVCAR8-ADR"]  (column name of log fold changes in LFC_table.csv)
   - annotation_table.csv
     The annotation_table has the following columns:
     Column1	rowname	                    gene_1	        gene_2	        gene_3	        guide_1	    guide_2	    guide_3
     0	        AAGCGAGT;AAGCGAGT;AAGCGAGT	dummyguide_1	dummyguide_1	dummyguide_1	AAGCGAGT	AAGCGAGT	AAGCGAGT

  To get the annotation table several steps are necessary:
  a) Download 'GSE154112_Three-way_ovarian_cancer_gRNA_screen_count_per_million_reads.xlsx'
     from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154112
     Save the first 4 columns (Key, mU6 (P3), U6 (P2), H1 (P1)) in the file raw_names.csv with the following names:
     "Key", "gene_1", "gene_2", "gene_3"

  b) Use the barcode_list.csv file which was required for running CombiPIPE:
    AAGCGAGT,1,dummyguide_1
    CTCTAGGT,2,dummyguide_2
    ACGTGAGA,3,CDK4_2
    ACTACTCG,4,CDK4_3
    AGTCCACA,5,DNMT1_1
    AGTTACAG,6,DNMT1_2
    ATAGGTGG,7,EGFR_1
    ATCACGCA,8,EGFR_2
    ATTCAGCC,9,ERBB2_1
    ATTCGGTT,10,ERBB2_2
    CACTTGTT,11,FGF2_1
    CAGCTTAT,12,FGF2_2
    CCTTCTCT,13,HDAC1_1
    CGAACTAG,14,HDAC1_2
    CTATCCAC,15,IKBKB_1
    CTCCAGTA,16,IKBKB_2
    GACGTTAC,17,MAP2K1_1
    GACTCTTA,18,MAP2K1_2
    GATTCGTT,19,MTOR_1
    GCAATCAA,20,MTOR_2
    TTACACGG,21,PIK3C3_1
    TTCAGGTC,22,PIK3C3_2
    GTCAGACT,23,POLA1_1
    GTCCTAGT,24,POLA1_2
    TTCGGCAA,25,TGFB1_1
    TTCTACGA,26,TGFB1_2
    TCCTATCA,27,TOP1_1
    TCCTGTAC,28,TOP1_2
    TCGGAACT,29,TUBA1A_1
    TCGTACTT,30,TUBA1A_2
    TCTCCGAT,31,TYMS_1
    TCTGTAGG,32,TYMS_3
'''

# Load data frames
raw_names = pd.read_csv("C://Gemi//Input//raw_names.csv", sep=";", index_col=0)
barcode_list = pd.read_csv("C://Gemi//Input//barcode_list.csv", sep=",", header=None, names=["guide", "number", "gene"])
LFC_raw_counts = pd.read_csv("C://Gemi//Input//LFC_raw_counts.csv", sep="\t")

raw_names[["guide_1", "guide_2", "guide_3", "rowname"]] = ""                                    # add four new columns

#  Map the guides to the genes by using the barcode.csv file
for row in range(0, len(raw_names)):
    for bar_row, bar_column in barcode_list.iterrows():
        if raw_names.iloc[row, 0] == bar_column["gene"]:
            raw_names.iloc[row, 3] = bar_column["guide"]
        if raw_names.iloc[row, 1] == bar_column["gene"]:
            raw_names.iloc[row, 4] = bar_column["guide"]
        if raw_names.iloc[row, 2] == bar_column["gene"]:
            raw_names.iloc[row, 5] = bar_column["guide"]

raw_names["rowname"] = raw_names["guide_1"] + ";" + raw_names["guide_2"] + ";" + raw_names["guide_3"]  # fill rowname column
raw_names["Key"] = raw_names["Key"].apply(lambda x: x.replace(",", "_"))                # change separator in key column

df_new = LFC_raw_counts.merge(raw_names, left_on='key', right_on='Key', how='left')     # merge LFCs to raw names df

# Create and save LFC_table and annotation_table
LFC_table = pd.DataFrame(df_new, columns=["rowname", "Cell_D26"])
annotation_table = pd.DataFrame(df_new, columns=["rowname", "gene_1", "gene_2", "gene_3", "guide_1", "guide_2", "guide_3"])
LFC_table.to_csv("C://Gemi//Input//LFC_table.csv", sep="\t", index=False)
annotation_table.to_csv("C://Gemi//Input//annotation_table.csv", sep="\t")

'''
    Tips for the annotation_table:
    - every gene should appear only one time, i.e. CDK1_1 and CDK1_2 should both be renamed to CDK1
    - no underscores are allowed in the gene names, e.g. dummyguide_1 needs to be renamed to dummyguide1

10. Run Gemi3 (see an example in run_example.py)
'''