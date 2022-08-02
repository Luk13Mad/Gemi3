# Gemi3

`Gemi3` is a Python package with a Rust core used for the analysis of 3D combinatorial CRISPR screens. </br>
The package provides out-of-the-box analysis solutions and graphics to illustrate the results.</br>
The Python interface allows integration into existing processing pipelines and also offers a high level of customisability.</br> 
The Rust core enables extremely fast execution and scalability with HPC clusters.</br> 
**A release to PyPI is planned.**

# Usage
Please refer to the example in the `example_run` directory.</br>
Example input data is also provided.

# Installation

1. Clone this repository</br>
2. Change directory to the downloaded repo</br>
3. Run the following command</br>
```
python setup.py install
```

# Citation
To cite `Gemi3` in your publications, please use: </br>

```
Bayesian Analysis for 3D combinatorial CRISPR screens
Lukas Madenach, Caroline Lohoff
bioRxiv 2022.06.02.494493; doi: https://doi.org/10.1101/2022.06.02.494493 
```

You can also use the following BibTeX entry:</br>

```
@article {Madenach2022.06.02.494493,
author = {Madenach, Lukas and Lohoff, Caroline},
title = {Bayesian Analysis for 3D combinatorial CRISPR screens},
elocation-id = {2022.06.02.494493},
year = {2022},
doi = {10.1101/2022.06.02.494493},
publisher = {Cold Spring Harbor Laboratory},
abstract = {Combinatorial CRISPR screens are a well-established tool for the investigation of genetic interactions in a high-throughput fashion. Currently advancements from 2D combinatorial CRISPR screens towards 3D combinatorial screens are made, but at the same time an easy-to-use computational method for the analysis of 3D combinatorial screens is missing. Here we propose a Bayesian analysis method for 3D CRISPR screens based on a well-established 2D CRISPR screen analysis protocol. Hoping to provide researchers with an out-of-the-box analysis solution, avoiding the need for time and resource intensive development of custom analysis protocols.Competing Interest StatementThe authors have declared no competing interest.},
URL = {https://www.biorxiv.org/content/early/2022/06/02/2022.06.02.494493},
eprint = {https://www.biorxiv.org/content/early/2022/06/02/2022.06.02.494493.full.pdf},
journal = {bioRxiv}
}

```

# Contact

Write a mail to `lukas.madenach{at}kitz-heidelberg.de` or just raise an issue.
