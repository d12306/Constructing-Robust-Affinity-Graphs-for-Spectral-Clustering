function clusts = SPClustering(A, num_clust)
clusts = zeros(size(A, 1), 1);
Clusters_gcut = gcut(A, num_clust);
for i = 1 : length(Clusters_gcut)
    clusts(Clusters_gcut{i}) = i;
end    

