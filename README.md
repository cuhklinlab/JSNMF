# JSNMF
Source codes and a demo of JSNMF are provided in this repository

1. JSNMF includes the main functions below:

  jsnmf_bpp_del_fs.mF: jsnmf implementation, a joint single-cell multi-omics integration approach based on semi-orthogonal nonnegative matrix factorization to dissect cellular heterogeneity

  preprocessing.R: data preprocessing script before running JSNMF.

  Wtrim.m: conduct Knn on cell-cell similairity matrix obtained from JSNMF.

  RAGI.m: evaluated the performance of the clustering methods using the Residual Average Gini Index (RAGI).
  
  demo_kidney.m: an example to use jsnmf method on mouse kidney data. It contains some downstream analysis tasks, such as clustering, visualization and so on.

2. Datasets and Examples

  Please refer to the vigenette with several examples for a quick guide to JSNMF.

3. Reference
