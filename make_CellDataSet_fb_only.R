#Make separate Cell Data Sets for forebrain and midbrain samples


##
## Merging the data
##
allCounts_fb <-  merge(fb_wt_1,fb_wt_2,by.x=0,by.y=0,all=TRUE)
allCounts_fb <-  merge(allCounts_fb,fb_ivs_1,by.x=1,by.y=0,all=TRUE)
allCounts_fb <-  merge(allCounts_fb,fb_ivs_2,by.x=1,by.y=0,all=TRUE)


##
## Writing the counts to file
##
rownames(allCounts_fb) <-  allCounts_fb[,1]
allCounts_fb <-  allCounts_fb[,-1]
write.csv(allCounts_fb,file='allCounts_fb.csv')

##
## Making the meta data
##
sData_fb <-  data.frame(
  'sampleName' = c('fb_wt_1','fb_wt_2','fb_ivs_1','fb_ivs_2'),
  'genotype' = c('wt','wt','ivs','ivs'),
  'region' = c(rep('fb',4)),
  'div' = rep(60,4),
  'genoReg' = c(rep('fbWt',2),(rep('fbIvs',2)))
)

allpData_fb <-  data.frame(
  'sampleName' = sub("-.*","",colnames(allCounts_fb)),
  'cellNum' = colnames(allCounts_fb)
)
allpData_fb <-  merge(allpData_fb,sData_fb,by.x=1,by.y=1,all.x=T,all.y=T,sort=F)
rownames(allpData_fb) = allpData_fb[,'cellNum']

## randomly sample the pData dataset
## to make sure that all the metadata looks good
## repeat this until you are convinced
allpData_fb[sample(1:dim(allpData_fb)[1],1),]


##
## filtering the data
library(Matrix)
## replace the NA values from the merge with 0
allCounts_fb[is.na(allCounts_fb)] = 0

## filtering out cells that do not have
## at least 100 reads.
validCells_fb <-  which(colSums(allCounts_fb) > 100)
length(validCells_fb)
##6397 cells total.
filteredCounts_fb <-  as(as.matrix(allCounts_fb[,names(validCells_fb)]), 'sparseMatrix')

## filtering out cells that do not fall within 2 sd of the ----------problems here with colSums
## mean number of reads
total_mRNAs_fb <-  Matrix::colSums(filteredCounts_fb)
upper_bound_fb <- 10^(mean(log10(total_mRNAs_fb))+
                     2*sd(log10(total_mRNAs_fb)))
lower_bound_fb <- 10^(mean(log10(total_mRNAs_fb)) -
                     2*sd(log10(total_mRNAs_fb)))

################ keeping only the cells in the upper and lower bound
filteredCounts_fb <-  filteredCounts_fb[,names(which(total_mRNAs_fb > lower_bound_fb & total_mRNAs_fb < upper_bound_fb))]

## filtering the genes -----------------------continue testing here
## requiring at least 4 read in at least
## 10 cells.
genesfData_fb <-  data.frame(rep("protein coding",dim(filteredCounts_fb)[1]),
                        row.names=rownames(filteredCounts_fb))
genesfData_fb[,'gene_short_name'] = rownames(filteredCounts_fb)

numCellExpressed_fb <-  Matrix::rowSums(filteredCounts_fb > 3)
genesfData_fb[,'numExpressed'] <-  numCellExpressed_fb
cellsPerGene_fb <-  Matrix::rowSums(filteredCounts_fb > 3)
expressedGenes_fb <-  names(which(cellsPerGene_fb > 10))
length(expressedGenes_fb)
#2737 genes expressed

##
## Loading into Monocle datframe
##
library(monocle)

##check that dimensions of the 3 input files are correct
dim(filteredCounts_fb[expressedGenes_fb,])
## 2737 rows 6103 columns
dim(allpData_fb[colnames(filteredCounts_fb),])
#should have tye same number of rows as Counts has.
#columns = 6103
dim(genesfData_fb[expressedGenes_fb,])
#should have the same number of rows as Counts has.
#rows =2737

HSMM_fb_1 <-  newCellDataSet(
  filteredCounts_fb[expressedGenes_fb,],
  phenoData = new("AnnotatedDataFrame",data=allpData_fb[colnames(filteredCounts_fb),]),
  featureData = new("AnnotatedDataFrame",data=genesfData_fb[expressedGenes_fb,]),
  expressionFamily=negbinomial())

HSMM_fb_1 <- detectGenes(HSMM_fb_1, min_expr = 4)
#sets global detection threshold to be used with this CellDataSet
HSMM_fb_1 <- estimateSizeFactors(HSMM_fb_1)
HSMM_fb_1 <- estimateDispersions(HSMM_fb_1)


qplot(total_mRNAs_fb[rownames(pData(HSMM_fb_1))], data=pData(HSMM_fb_1), color=genotype,geom="density")+
  geom_vline(xintercept=lower_bound_fb) +
  geom_vline(xintercept=upper_bound_fb)

qplot(total_mRNAs_fb[rownames(pData(HSMM_fb_1))], data=pData(HSMM_fb_1), color=sampleName,geom="density")+
  geom_vline(xintercept=lower_bound_fb) +
  geom_vline(xintercept=upper_bound_fb)




