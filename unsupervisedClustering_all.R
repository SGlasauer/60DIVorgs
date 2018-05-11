#unsupervised clustering of cells - all 60DIV samples

#make table of variance of a gene´s expression across cells
disp_table_all <- dispersionTable(HSMM1)
#filter out cells with low mean expression
unsup_clustering_genes_all <- subset(disp_table_all, mean_expression >= 0.1)

#setOrderingFilter to marks genes that will be used for clustering
HSMM1 <- setOrderingFilter(HSMM1, unsup_clustering_genes_all$gene_id)

#plot dispersion against average expression of a gene
plot_ordering_genes(HSMM1)


plot_pc_variance_explained(HSMM1, return_all = F) # norm_method='log'
#this gives me a warning: did not converge--results might be invlaid!
#Dave Tangs blog had the same error, so I guess its fine, but keep this in mind!


#now reduce dimensions. I chose 5 based in the plot explaining variance.
HSMM1 <- reduceDimension(HSMM1, max_components = 2, num_dim = 5,
                        reduction_method = 'tSNE', verbose = T)
#this gave me a Warning message:
##In doTryCatch(return(expr), name, parentenv, handler) :
##restarting interrupted promise evaluation

#to cluster the cels: requires how many clusters you want as input
HSMM1 <- clusterCells(HSMM1, num_clusters = 15)
##Message: Distance cutoff calculated to 4.93999

plot_cell_clusters(HSMM1)

##check whether previously assigned cell types cluster together
plot_cell_clusters(HSMM1, color = "CellType")
plot_cell_clusters(HSMM1, color = "genoReg")
plot_cell_clusters(HSMM1, color = "region")
plot_cell_clusters(HSMM1, color = "sampleName")
