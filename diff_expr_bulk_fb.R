#Differential gene expression on whole forebrain dataset


##you can first select the genes that you want to do differential analysis on
##but here I am going to use all the genes in the forebrain dataset

#check how many genes the dataset has
length(row.names(fData(HSMM_fb_1)))
##2737 genes will be tested

diff_test_res_fb <- differentialGeneTest(HSMM_fb_1,
                                      fullModelFormulaStr = "~genotype")

# Select genes that are significant at an FDR < 10%
sig_genes_fb <- subset(diff_test_res_fb, qval < 0.1)

nrow(sig_genes_fb)
##1537 genes are differentially expressed between ctrl and ivs
##that seems like a lot

sig_genes_fb_ordered <- dplyr::arrange(sig_genes_fb, pval, qval)

write.csv(sig_genes_fb_ordered, file='sig_genes_fb_ordered.csv')

##plot a few selected genes that appear on top of the list

head(sig_genes_fb_ordered)

diff_genes_selection0 <- HSMM_fb_1[row.names(subset(fData(HSMM_fb_1),
                                                    gene_short_name %in% c("MEG3", "MTRNR2L1", "TMSB4X"))),]
plot_genes_jitter(diff_genes_selection0, grouping = "genotype", ncol= 3)
plot_genes_jitter(diff_genes_selection0, grouping = "sampleName", ncol= 3)


diff_genes_selection1 <- HSMM_fb_1[row.names(subset(fData(HSMM_fb_1),
                                      gene_short_name %in% c("MAP2", "MAPT"))),]
plot_genes_jitter(diff_genes_selection1, grouping = "genotype", ncol= 2)
plot_genes_jitter(diff_genes_selection1, grouping = "sampleName", ncol= 2)

diff_genes_selection2 <- HSMM_fb_1[row.names(subset(fData(HSMM_fb_1),
                                                    gene_short_name %in% c("PAX6", "NES"))),]
plot_genes_jitter(diff_genes_selection2, grouping = "genotype", ncol= 2)
plot_genes_jitter(diff_genes_selection2, grouping = "sampleName", ncol= 2)

diff_genes_selection3 <- HSMM_fb_1[row.names(subset(fData(HSMM_fb_1),
                                                    gene_short_name %in% c("SFRP1", "SFRP2"))),]
plot_genes_jitter(diff_genes_selection3, grouping = "genotype", ncol= 2)
plot_genes_jitter(diff_genes_selection3, grouping = "sampleName", ncol= 2)

diff_genes_selection4 <- HSMM_fb_1[row.names(subset(fData(HSMM_fb_1),
                                                    gene_short_name %in% c("KIF3A", "KIF5C"))),]
plot_genes_jitter(diff_genes_selection4, grouping = "genotype", ncol= 2)
plot_genes_jitter(diff_genes_selection4, grouping = "sampleName", ncol= 2)

diff_genes_selection5 <- HSMM_fb_1[row.names(subset(fData(HSMM_fb_1),
                                                    gene_short_name %in% c("TIMP3", "SOX11"))),]
plot_genes_jitter(diff_genes_selection5, grouping = "genotype", ncol= 2)
plot_genes_jitter(diff_genes_selection5, grouping = "sampleName", ncol= 2)


####try diferential gene expression test between clusters -with a few genes first.
## I cannot find a "Cluster" column in the pData, therefore I am not sure whether it will work.


fb_genes_to_test <- row.names(subset(fData(HSMM_fb_1),
                                 gene_short_name %in% c("PAX6", "MAP2", "KIF3A",
                                                        "KIF5C", "MAPT","SFRP1",
                                                        "SFRP2",  "GFAP",  "NCAM",
                                                        "SOX4",  "SOX11", "NES",
                                                        "FABP7", "FOXP1", "PCP4",
                                                        "OTX2",  "TCF4", "RGS2",
                                                        "TIMP3", "JUN", "SLC1A3")))



diff_test_res_fb_clusters_selecetedMarkers <- differentialGeneTest(HSMM_fb_1,
                                         fullModelFormulaStr = "~Cluster")

# Select genes that are significant at an FDR < 10%
sig_genes_fbclusters_selectedMarkers <- subset(diff_test_res_fb, qval < 0.1)

nrow(sig_genes_fbclusters_selectedMarkers)
##1537 genes are differentially expressed between ctrl and ivs

