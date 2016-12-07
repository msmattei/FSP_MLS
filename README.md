# FSP_MLS
First step project (MLS master program) with Marc Robinson-Rechavi.

data_extraction.R: R code to extract data information from biomaRt and Bgee. Combine them into a merged database (saved as tissue data) with expression data (RPKM coming from Bgee), gene id (from biomaRt) and protein domain information (from Erich Bornberg).

tissue_specificity_score.R: R code computing the tissue specificity score for each experiment
