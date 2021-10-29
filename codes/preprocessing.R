# preprocess mouse kidney data
library(data.table)
library(dplyr)
library(Seurat)
# setwd("data/")
# preprocessing kidney and select hvg features using Seurat package
# rna.csv, atac.csv: scRNA count data and scATAC data
# barcodes.csv: barcodes of cells
# genes.csv: genes for rna data
# peaks.csv: regions for atac data
# these data above are included in the "Kidney_sciCAR_data.mat" file, and need to save .csv files, respectively

rna<-read.table("rna.csv",sep = ',',header = FALSE)
barcodes<-read.csv("barcodes.csv",sep = ',',header = FALSE)
genes<-read.csv("genes.csv",sep = ',',header = FALSE)
row.names(rna) <- t(genes)
names(rna) =t(barcodes)

atac<-read.table("atac.csv",sep = ',',header = FALSE)
peaks<-read.csv("peaks.csv",sep = ',',header = FALSE)
row.names(atac) <- t(peaks)
names(atac) =t(barcodes)

kidney <- CreateSeuratObject(counts = rna, project = "kidney", min.cells = 10, min.features = 500)
kidney <- NormalizeData(kidney, normalization.method = "LogNormalize", scale.factor = 10000)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 5000)
top5k <- head(VariableFeatures(kidney), 5000)
write.csv(kidney@assays[["RNA"]]@data[top5k,],"rna_norm5k.csv")

memory.limit(size=56000) 
kidney <- CreateSeuratObject(counts = atac, project = "kidney", min.cells = 5, min.features = 200)
kidney <- NormalizeData(kidney, normalization.method = "LogNormalize", scale.factor = 10000)
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 10000)
top10k <- head(VariableFeatures(kidney), 10000)
write.csv(kidney@assays[["kidney"]]@data[top10k,],"atac_norm10k.csv")



