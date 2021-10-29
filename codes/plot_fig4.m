% prepare for Figure 4 on mouse skin data---an illustrative example

% marker_list.csv was downloaded from the original publication of "SHARE-seq", and used to extract marker genes
marker_list = importdata('data/mouse_skin/marker_list.csv');
mark_value = marker_list.data;
celltype = marker_list.textdata(1,2:end);
mark_genes = marker_list.textdata(2:end,1);
clear marker_list
% import W1 obtained from JSNMF
load('data/mouse_skin/skin_W1.mat'); 
% import features (genes) 
load('data/mouse_skin/skin_5k_10k.mat','genes'); 

[C, ind] = sort(mark_value,'descend');
%index 16 corresponds to 'Dermal Fibrobalst' cell type, and index 9: 'Hair Shaft-Cuticle/Cortex  '
index = [16,9];  
for i = 1:length(index)
    tmp = find(mark_value(:,index(i)) == 1);
    fc_genes{i} = mark_genes(tmp);
end

fc_celltype = celltype(index);
% In W1, the column 3 and 11 correspond to 'Dermal Fibrobalst' and 'Hair Shaft-Cuticle/Cortex' cell types
idx_in_w1 = [3,11];
top_per = [];
featureRankingPlot_1(W1,genes,fc_genes,[],[],idx_in_w1,top_per);



