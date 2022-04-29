############################################### Example usage of Gemi3 #################################################

import Gemi3.g3datahandling as g3d
import Gemi3.Model1 as md1
import Gemi3.Score as sc

# Load data
#LFC table can be generated with g3datahandling.calc_observed_LFC()
LFC_table = "~/Desktop/rust_projects/Gemi3/example_run/LFC_table_notnan_Zhou.csv"
annotation_table = "~/Desktop/rust_projects/Gemi3/example_run/annotation_table_Zhou.csv"
samples = ["OVCAR8-ADR"]  # column name with log2 fold changes in LFC_table

# Run g3datahandling
mydata = g3d.G3D(LFC_table,annotation_table,samples,True)

# Initialize statistical model
'''
The goal is to reduce the mean absolute error after every iteration. 
If the data contains only one cell line, prio_tau_x and prio_tau_y should be bigger (and therefor SD smaller) 
than prio_tau_r and prio_tau_s,because we expect the sample-independent effects (x) to be small.
It may take a while to find the right initial parameters. 
A tip is to look at the values of mymodel after the convergence failed. Setting the parameters to the actual values may
help to get the right dimensions and differences between the parameters (example: prio_tau_x = 1/np.std(mymodel.x)).
'''
mymodel = md1.Model1(mydata,["dummyguide1","dummyguide2"],verbose=True,force=False,
                     prio_tau_x=5.0, prio_tau_xx=5.0,
                     prio_tau_r=10.0, prio_tau_rr=10.0,
                     prio_tau_y=10.0, prio_tau_yy=10.0,
                     prio_tau_s=10.0, prio_tau_ss=10.0,
                     prio_mu_x=1.0, prio_mu_xx=1.0,
                     prio_mu_y=0.0, prio_mu_yy=0.0,
                     prio_mu_r=1.0, prio_mu_rr=1.0,
                     prio_mu_s=0.0, prio_mu_ss=0.0
                     )

# Prime model by initializing x,y,r,s and tau values and setting the hash-tables
mymodel.prime()

# Inference using CAVI implemented in rust
mymodel.start_inference_rust()

# Compute scores
'''
- use the class 'Score' to evaluate genetic interaction
- define at least one positive control gene (to remove interactions involving individually lethal genes)
- define at least one negative control gene (to calculate p-values and FDRs) 
- all data frames and plots are automatically saved in the folder 'Output'
'''

scores = sc.Score(model=mymodel, nc_genes=['CDK4'], pc_genes="EGFR",
                  pc_weight=0.5, verbose=True,plotting=True)


