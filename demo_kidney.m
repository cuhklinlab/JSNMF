%% run JSNMF algorithm on mouse kidney dataset
% load kindey data, including single-cell transciptome and epigenomic proiles
load('data/kidney_5k_10k.mat')
load('data/housekeeping_genes.mat')
load('data/Kidney_sciCAR_data.mat','RNA')
% data preprocessing using R package "Seurat 4.0", preprocessed data was used in the analysis below
% we provide script "preprocessing.R" to conduct data preprocess, normalization and hvg selection steps in R enviroment

% preparing for computing RAGI using ind1 ind2 data and so on 
[~,ix1,ix2]=intersect(share_bd,RNA.Cells,'stable');
genes = RNA.Features; data = RNA.data;
data = data(:,ix2);
sM = sum(data,2); zero_row = find(sM==0);
data(zero_row,:) = []; genes(zero_row) = [];
data = data./repmat(sum(data),size(data,1),1)*10000;
data = log(data+1);

label_name = label;
num_clu = length(unique(label_name));
tag = unique(label_name); tag = cellstr(tag);
[~, label] = ismember(label_name,tag);

% marker genes are given by the original publication; house-keeping genes are given in the website "http://www.housekeeping.unicamp.br/"
marker_genes = ["Slc12a3";'Trpm6';'Abca13';'Klhl3';'Wnk1';'Tsc22d1';'Cadps2';...
    'Egfem1';'Calb1';'Kl';'Temem72';'Dach1';'Ptprm';'Sgms2';'Tox3';'Frmpd4';...
    'Rbms3';'Fxyd4';'Dnm3';'Cacnb2';'Pde1c';'Pde8b'];
[~,~,ind1] = intersect(marker_genes,genes,'stable');
[~,~,ind2] = intersect(hp_genes,genes,'stable');
data = sparse(data);
clear ix1 ix2 sM zero_row tag share_bd RNA peaks marker_genes hp_genes genes

% parameter seletion 
[alpha,gamma,Inits] = parameter_selection(X1,X2,label);
% start timing
tic
[W1,W2,H1,H2,S,iter,objs] = jsnmf(X1,X2,alpha,gamma,Inits);
disp('JSNMF runtime:');
toc

A = Wtrim(S,50); % KNN on S, k is insensitive to the performance, here we suggested k=50
[clust,~,~] = getNCluster(A,num_clu,0,3,20); % implement louvain clustering
if length(unique(clust))== num_clu
   [ac, nmi_value, ~] = CalcMetrics(label, clust);
   [ave_mk_gini, ave_hk_gini, difgini] = RAGI(data,ind1,ind2,clust);
else
% using spectral clustering instead if number of clusering doesn't equals
% to the pointed number; in general we can obtain the pointed number by adjusting the resolution parameter
   [~, clust, ~] = SpectralClustering(A, num_clu);
   [ac, nmi_value, ~] = CalcMetrics(label, clust);   
   [ave_mk_gini, ave_hk_gini, difgini] = RAGI(data,ind1,ind2,clust);
end

% visualization on S
addpath('umapFileExchange/umap') % load UMAP package
addpath('umapFileExchange/util')

term = '.'; 
colors = generateColors(max(length(unique(label_name))),length(unique(term)));

D_jsnfm = 1-S; D_jsnfm = D_jsnfm-diag(diag(D_jsnfm));
[reduction_jsnmf, ~, ~, ~] = run_umap(D_jsnfm,'metric','precomputed'); 
gscatter(reduction_jsnmf(:,1),reduction_jsnmf(:,2),label_name,colors,[],4); 
set(gca,'xtick',[],'ytick',[]);
title('UMAP for JSNMF')

legend('Location','westoutside','Box','off','FontSize',9.5); 
legendmarkeradjust(16)
legend('boxoff') 


