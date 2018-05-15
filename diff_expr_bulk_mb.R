#Differential gene expression on whole midbrain dataset


##you can first select the genes that you want to do differential analysis on
##but here I am going to use all the genes in the midbrain dataset

#check how many genes the dataset has
length(row.names(fData(HSMM_mb_1)))
##2801 genes will be tested

diff_test_res_mb <- differentialGeneTest(HSMM_mb_1,
                                         fullModelFormulaStr = "~genotype")

# Select genes that are significant at an FDR < 10%
sig_genes_mb <- subset(diff_test_res_mb, qval < 0.1)

nrow(sig_genes_mb)
##1710 genes are differentially expressed between ctrl and ivs
##that seems like a lot

sig_genes_mb_ordered <- dplyr::arrange(sig_genes_mb, pval, qval)

write.csv(sig_genes_mb_ordered, file='sig_genes_mb_ordered.csv')


##pot a few selected genes that appear on top of the list

head(sig_genes_mb_ordered)

diff_genes_mb_selection0 <- HSMM_mb_1[row.names(subset(fData(HSMM_mb_1),
                                                    gene_short_name %in% c("MEG3", "NNAT", "TMSB10"))),]
plot_genes_jitter(diff_genes_mb_selection0, grouping = "genotype", ncol= 3)
plot_genes_jitter(diff_genes_mb_selection0, grouping = "sampleName", ncol= 3)


diff_genes_mb_selection1 <- HSMM_mb_1[row.names(subset(fData(HSMM_mb_1),
                                                       gene_short_name %in% c("MAP2", "MAPT"))),]
plot_genes_jitter(diff_genes_mb_selection1, grouping = "genotype", ncol= 2)
plot_genes_jitter(diff_genes_mb_selection1, grouping = "sampleName", ncol= 2)

