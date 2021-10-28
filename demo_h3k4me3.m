%% run JSNMF algorithm on mouse brain (RNA+H3K4me3, paired-TAG) data
load('data/h3k4me3_5k_10k.mat')
num_clu = length(unique(label));
marker_genes = ["Dnah12"; "Cfap299";'Flt1'; 'Slco1a4'; 'Inpp5d'; 'Hexb'; 'Slc1a3'; 'Atp1a2'; 'Slc1a2'; 'Gpc5'; ...
    'Prr5l'; 'Plp1'; 'Rnf220'; 'Lhfpl3'; 'Cacng4'; 'Erbb4'; 'Kcnc2'; 'Reln'; 'Grin3a'; 'Adarb2'; 'Fam19a2'; ...
    'Ahcyl2'; 'Gm32647'; 'Shisa6'; 'Epha6'; 'Galnt14'; 'Hs3st4'; 'Tshz2'; 'Fam19a1'; 'Il1rapl2'; 'Cdh12'; 'Rgs6';...
    'Lingo2'; 'Prr16'; 'Zfp804b'; 'Pex5l'; 'Gm26883'; 'Cdh18'; 'Gm28928'; 'Galntl6'; 'Spag16'; 'Cfap43'; 'Wdr49';...
    'Ebf1'; 'Mecom'; 'Ptprb'; 'Tgfbr1'; 'Zfhx3'; 'Apbb1ip'; 'Gpc5'; 'Plpp3'; 'Slc1a2'; 'Igsf8'; 'St18'; 'Mog'; 'Mag';...
    'St18'; 'Mbp'; 'Vcan'; 'Tnr'; '6030443J06Rik'; 'Nxph1'; '6330411D24Rik'; 'Kcnmb2'; 'Grip1'; 'Gm45341'; 'Gm45455';...
    'Grip1'; 'Trpm3'; 'Rfx3'; 'Dgkh'; 'Hs3st4'; 'Ryr3'; 'Hs6st3'; 'Gm2164'; 'Grik4'; 'Cpne7'; 'Gm2164'; 'Ryr3';...
    'Hs6st3'; 'Sdk1'; 'Foxp2'; 'Garnl3hrm2'; 'Grm8'; 'Vwc2l'; 'Olfm3'; 'Grik1 Gm2164'; 'Gm28928'; 'Sgcz'; 'Prr16'; ...
    'Chrm3'; 'Pdzrn3'; 'Kcnq5'; 'Unc5d'; 'Car10'; 'Pcdh15'; 'Nrg1']; 
marker_genes = unique(marker_genes);
load('data/housekeeping_genes.mat')
housekeeping_genes = hp_genes;

load('data/h3k4me3_rna_RAGI.mat')
genes = genes_h3k4me3; data = rna_h3k4me3;
[~,~,ind1] = intersect(marker_genes,genes,'stable');
[~,~,ind2] = intersect(housekeeping_genes,genes,'stable');

% parameter seletion 
[alpha,gamma,Inits] = parameter_selection(X1,X2,label);
% start timing
tic
[W1,W2,H1,H2,S,iter,objs] = jsnmf(X1,X2,alpha,gamma,Inits);
disp('JSNMF runtime:');
toc

A = Wtrim(S,50); % KNN on S
[clust,~,~] = getNCluster(A,num_clu,0,3,20); % implement louvain clustering
if length(unique(clust))== num_clu
   [ac, nmi_value, ~] = CalcMetrics(label, clust);
   [ave_mk_gini, ave_hk_gini, difgini] = RAGI(data,ind1,ind2,clust);
else
   [~, clust, ~] = SpectralClustering(A, num_clu); % using spectral clustering instead if number of clusering doesn't equals to the pointed number
   [ac, nmi_value, ~] = CalcMetrics(label, clust);   
   [ave_mk_gini, ave_hk_gini, difgini] = RAGI(data,ind1,ind2,clust);
end

% visualization on S
addpath('umapFileExchange/umap') % load UMAP package
addpath('umapFileExchange/util')

term = '.'; 
label_name = label_h3k4me3(:,2);
colors = generateColors(max(length(unique(label_name))),length(unique(term)));

D_jsnfm = 1-S; D_jsnfm = D_jsnfm-diag(diag(D_jsnfm));
[reduction_jsnmf, ~, ~, ~]=run_umap(D_jsnfm,'metric','precomputed','min_dist',0.68,'n_neighbors',12); %,'min_dist',0.68,'n_neighbors',12
gscatter(reduction_jsnmf(:,1),reduction_jsnmf(:,2),label_name,colors,[],4); 
set(gca,'xtick',[],'ytick',[]);
title('UMAP for JSNMF')

legend('Location','westoutside','Box','off','FontSize',9.5); 
legendmarkeradjust(16)
legend('boxoff') 


