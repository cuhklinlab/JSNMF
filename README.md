# JSNMF
Source codes and a demo of JSNMF are provided in this repository

1. JSNMF includes the main functions below:

  jsnmf.mF: jsnmf implementation, a joint single-cell multi-omics integration approach based on semi-orthogonal nonnegative matrix factorization to dissect cellular heterogeneity

  preprocessing.R: data preprocessing script before running JSNMF.
  
  parameter_selection.m: select parameter of JSNMF(i.e. two regularization parameters alpha and gamma).
  
  SNF.m: similarity network fusion script provied in the original publication "Similarity network fusion for aggregating data types on a genomic scale". 
  
  Wtrim.m: conduct Knn on cell-cell similairity matrix obtained from JSNMF.

  RAGI.m: evaluated the performance of the clustering methods using the Residual Average Gini Index (RAGI).
  
  plot_fig4.m: plot figure 4 in the body text
  
  demo_kidney.m: an example to use jsnmf method on mouse kidney data (sci-CAR). It contains some downstream analysis tasks, such as clustering, visualization and so on.
  
  demo_h3k4me3.m: an example to use jsnmf method on mouse brain data (scRNA + H3K4me3, Paired-Tag). It contains some downstream analysis tasks, such as clustering, visualization and so on.

2. Datasets and Examples

  Please refer to two examples(jsnmf_mouse brain_H3K4me3_paired-Tag.mlx and jsnmf_mouse kidney_sci-CAR.mlx) for a quick guide to JSNMF.

3. Reference
