#unsupervised clustering of cells - only midbrain 60DIV samples

#make table of variance of a gene´s expression across cells
disp_table_mb <- dispersionTable(HSMM_mb_1)
#filter out cells with low mean expression
unsup_clustering_genes_mb <- subset(disp_table_mb, mean_expression >= 0.1)

#setOrderingFilter to marks genes that will be used for clustering
HSMM_mb_1 <- setOrderingFilter(HSMM_mb_1, unsup_clustering_genes_mb$gene_id)

#plot dispersion against average expression of a gene
plot_ordering_genes(HSMM_mb_1)


plot_pc_variance_explained(HSMM_mb_1, return_all = F) # norm_method='log'
#this gives me a warning: did not converge--results might be invlaid!
#Dave Tangs blog had the same error, so I guess its fine, but keep this in mind!


#now reduce dimensions. I chose 6 based in the plot explaining variance.
HSMM_mb_1 <- reduceDimension(HSMM_mb_1, max_components = 2, num_dim = 6,
                             reduction_method = 'tSNE', verbose = T)


#to cluster the cels: requires how many clusters you want as input
HSMM_mb_1 <- clusterCells(HSMM_mb_1, num_clusters = 10)
##Message: Distance cutoff calculated to 4.476947

plot_cell_clusters(HSMM_mb_1)

plot_cell_clusters(HSMM_mb_1, color = "genotype")
plot_cell_clusters(HSMM_mb_1, color = "sampleName")


