# Pathonoia

This repository contains the software (Python Package and Jupyter Notebook) descibed in:

- Liebhoff, AM. et al. **Pathogen detection in RNA-seq data with Pathonoia**. *BMC Bioinformatics* 24, 53 (2023). https://doi.org/10.1186/s12859-023-05144-z

The `pathonoia_notebook` Jupyter Notebook in the top directory shows the use of the Python package for a whole dataset.
To reproduce the results published in the article and getting Pathonoia to work with minimal efforts, we provide a precompiled Kraken database here:

[oasis.ims.bio/kraken2_k31.zip](https://oasis.ims.bio/kraken2_k31.zip) (March 2019)

An updated version of the index is provided here:

[oasis.ims.bio/kraken2_k31_2205.zip](https://oasis.ims.bio/kraken2_k31_2205.zip) (May 2022)

The Jupyter Notebook also explains the setup of Kraken index for Pathonoia if an adjusted version is desired. 
Be aware that Pathonoia is build for Kraken indexes that have k=l, i.e. minimizer size = k-mer size. 
Another index may work, but results and specificity cannot be guaranteed.

The folder `downstream_analysis` contains the Jupyter Notebooks with the code for the results presented in the paper and additionally a template for the analysis of a dataset analysed by Pathonoia (especially, after using the pathonoia_notebook in the main directory).

The folder `dataset` contains simulated sequences for an example dataset. The full dataset was also used for the publication. The data here, is a subset (files reduced in size, fewer sequences) to speed up the example and reduce the size of this repository. We thank Simon H.Ye for providing the simulated samples:
- Simon, H. Y., Siddle, K. J., Park, D. J., & Sabeti, P. C. (2019). Benchmarking metagenomics tools for taxonomic classification. *Cell*, 178(4), 779-794.
