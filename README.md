# JSNMF enables effective and accurate integrative analysis of single-cell multiomics data

Source codes and a demo of JSNMF are provided in this repository.

Running environment：``MATLAB R2019b`` or later. 

The external functions used in the manuscript and JSNMF can be found in the folder ``external/``.

# 1. JSNMF includes the main functions below:

  ``jsnmf.m``: jsnmf implementation, a joint single-cell multi-omics integration approach based on semi-orthogonal nonnegative matrix factorization to dissect cellular heterogeneity

  ``preprocessing.R``: data preprocessing script before running JSNMF.
  
  ``parameter_selection.m``: parameters selection for JSNMF(i.e. two hyperparameters alpha and gamma).
  
  ``SNF.m``: similarity network fusion script provided in the original publication "[Similarity network fusion for aggregating data types on a genomic scale](https://www.nature.com/articles/nmeth.2810)". 
  
  ``Wtrim.m``: conduct KNN on cell-cell similairity matrix obtained from JSNMF.

  ``RAGI.m``: evaluate the performance of the clustering methods using the Residual Average Gini Index (RAGI).
  
  ``plot_fig4.m``: plot figure 4 in the body text.
  
  ``sameH_jsnmf.m``: implement JSNMF (same H) on mouse kidney data.
  
  ``jsnmf_batch_correct_pari5.m``: an extension of JSNMF to correct batch effect on all five Paired-Tag datasets. In this setting, there are six modalities in total and the modality of RNA is shared in the five single-cell experiments. 
  
  ``jsnmf_3mod.m``: an extension of JSNMF to analyze single-cell multi-omics data with more than two modalities measured from the same single cell. A scNMT-seq dataset are used as an illustrative example.
  
  ``parameter_selection_batcheffect_pair5.m``, ``parameter_select_for_same_H.m`` and ``pa_sel.m``: parameter selection for ``jsnmf_batch_correct_pari5.m``, ``sameH_jsnmf.m`` and ``jsnmf_3mod.m``, respectively.
  
# 2. Datasets and Examples

  The datasets used in this manuscript have been uploaded.
  
  Please refer to two examples for a quick guide to JSNMF:

  ``jsnmf_mouse brain_H3K4me3_paired-Tag.mlx``: an example to use jsnmf method on mouse brain data (scRNA + H3K4me3, Paired-Tag). It contains some downstream analysis tasks, such as clustering, visualization and so on.

  ``jsnmf_mouse kidney_sci-CAR.mlx``: an example to use jsnmf method on mouse kidney data (sci-CAR). It contains some downstream analysis tasks, such as clustering, visualization and so on.

  We also include a tutorial of JSNMF (same H) in ``jsnmf_sameH.mlx``.
  
  JSNMF is generalizable to integrate multiple single-cell multi-omics experiments and single-cell experiments with more than two molecular modalities profiled in the same cell. In the first example, we extended JSNMF to integrate all five Paired-Tag datasets (RNA + H3K27me3, RNA + H3K4me3, RNA + H3K9me3, RNA + H3K27ac and RNA + H3K4me1) simultaneously. In the second example, we extended JSNMF to analyze single-cell multi-omics data with more than two modalities measured from the same single cell. We analyzed a scNMT-seq dataset generated from the developing mouse embryo, where gene expression, chromatin accessibility and methylation were profiled simultaneously in the same cell. 
  The tutorials for JSNMF to integrate multiple single-cell multi-omics experiments and more molecular modalities are presented in ``jsnmf_for_five_paired_tag.mlx`` and ``jsnmf_for_three_mod.mlx``, respectively.
  
  The tutorial for inferring cell-type-specific region-gene associations is presented in ``example/Gene_Assoc_Tutorial.md``.
  
  If you have any questions about the source code, please feel free to contact me: chonghua_1983@126.com. The python implementation of JSNMF are also provided in "https://github.com/cuhklinlab/JSNMF_py". 

# 3. Reference

Ma, Y., Sun, Z., Zeng, P., Zhang, W., & Lin, Z. (2022). [JSNMF enables effective and accurate integrative analysis of single-cell multiomics data](https://academic.oup.com/bib/article/23/3/bbac105/6563185). *Briefings in Bioinformatics*, 23(3), bbac105.