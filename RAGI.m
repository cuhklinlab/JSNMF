%% statistic the population size of different clusters to computer GINI indices
% since no cell type labels are provided in this dataset, e.g. 10X pbmc data
function [ave_marker_gini,ave_housekeeping_gini,diffgini] = RAGI(data,ind1,ind2,clust)
data = full(data);
num_clu = length(unique(clust));
%[~,~,ind1] = intersect(marker_genes,genes,'stable');
df_matrix_marker = data(ind1,:);
%[~,~,ind2] = intersect(housekeeping_genes,genes,'stable');
df_matrix_housekeeping = data(ind2,:);

% compute the mean of the expression values in all cells for each cluster
% to quantify the enrichment of each gene in each cluster of cells
df_marker = [clust, df_matrix_marker'];
pop = accumarray(df_marker(:,1),1);
ave_marker = zeros(num_clu,length(ind1));
marker_gini = zeros(length(ind1),1);

for i = 2:length(ind1)+1
    ave_marker(:,i) = accumarray(df_marker(:,1), df_marker(:,i)) ./ accumarray(df_marker(:,1),1);
    marker_gini(i-1) = gini(pop,ave_marker(:,i));
end

df_housekeeping = [clust, df_matrix_housekeeping'];
ave_housekeeping = zeros(num_clu,length(ind2));
housekeeping_gini = zeros(length(ind2),1);

for j = 2:length(ind2)+1
    ave_housekeeping(:,j) = accumarray(df_housekeeping(:,1), df_housekeeping(:,j)) ./ accumarray(df_housekeeping(:,1),1);
    housekeeping_gini(j-1) = gini(pop,ave_housekeeping(:,j));
end   

ave_marker_gini = mean(marker_gini);
ave_housekeeping_gini = mean(housekeeping_gini);
diffgini = ave_marker_gini-ave_housekeeping_gini;

disp(sprintf('ave_marker_gene: %0.4f\tave_hk_gini:%0.4f\tdiff_gini:%0.4f\t', ave_marker_gini, ave_housekeeping_gini, diffgini))




