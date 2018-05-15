#differential gene expression between fb and wt samples on whole dataset

#check how many genes the dataset has
length(row.names(fData(HSMM1)))
##3648 genes will be tested

diff_test_res_FBvsMB <- differentialGeneTest(HSMM1,
                                         fullModelFormulaStr = "~region")

# Select genes that are significant at an FDR < 10%
sig_genes_FBvsMB <- subset(diff_test_res_FBvsMB, qval < 0.1)

nrow(sig_genes_FBvsMB)
##2856 genes are differentially expressed between ctrl and ivs
##that seems like a lot

sig_genes_FBvsMB_ordered <- dplyr::arrange(sig_genes_FBvsMB, pval, qval)

write.csv(sig_genes_FBvsMB_ordered, file='sig_genes_FBvsMB_ordered.csv')


##plot a few selected genes that appear on top of the list

head(sig_genes_FBvsMB_ordered)
#top 5 genes: FABP7, FOXP1, LHX5-AS1, MALAT1, NR2F1, PCSK1

diff_genes_FBvsMB_selection0 <- HSMM1[row.names(subset(fData(HSMM1),
                                                  gene_short_name %in% c("FABP7", "FOXP1"))),]
plot_genes_jitter(diff_genes_FBvsMB_selection0, grouping = "region", ncol= 2)
plot_genes_jitter(diff_genes_FBvsMB_selection0, grouping = "sampleName", ncol= 1)
#####


diff_genes_FBvsMB_SULF2 <- HSMM1[row.names(subset(fData(HSMM1),
                                                  gene_short_name %in% ("SULF2"))),]
plot_genes_jitter(diff_genes_FBvsMB_SULF2, grouping = "region", ncol= 1)
plot_genes_jitter(diff_genes_FBvsMB_SULF2, grouping = "sampleName", ncol= 1)

plot_cell_clusters(HSMM1, color_by = "region", markers = "SULF2")

