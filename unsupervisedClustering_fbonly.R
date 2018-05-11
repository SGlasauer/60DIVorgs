#unsupervised clustering of cells - only forebrain 60DIV samples

#make table of variance of a gene´s expression across cells
disp_table_fb <- dispersionTable(HSMM_fb_1)
#filter out cells with low mean expression
unsup_clustering_genes_fb <- subset(disp_table_fb, mean_expression >= 0.1)

#setOrderingFilter to marks genes that will be used for clustering
HSMM_fb_1 <- setOrderingFilter(HSMM_fb_1, unsup_clustering_genes_fb$gene_id)

#plot dispersion against average expression of a gene
plot_ordering_genes(HSMM_fb_1)


plot_pc_variance_explained(HSMM_fb_1, return_all = F) # norm_method='log'
#this gives me a warning: did not converge--results might be invlaid!
#Dave Tangs blog had the same error, so I guess its fine, but keep this in mind!


#now reduce dimensions. I chose 8 based in the plot explaining variance.
HSMM_fb_1 <- reduceDimension(HSMM_fb_1, max_components = 2, num_dim = 8,
                        reduction_method = 'tSNE', verbose = T)
#this gave me a Warning message:
##In doTryCatch(return(expr), name, parentenv, handler) :
##restarting interrupted promise evaluation

#to cluster the cels: requires how many clusters you want as input
HSMM_fb_1 <- clusterCells(HSMM_fb_1, num_clusters = 15)
##Message: Distance cutoff calculated to 4.93999

plot_cell_clusters(HSMM_fb_1)

##check whether previously assigned cell types cluster together
plot_cell_clusters(HSMM_fb_1, color = "CellType")
plot_cell_clusters(HSMM_fb_1, color = "genotype")
plot_cell_clusters(HSMM_fb_1, color = "sampleName")


