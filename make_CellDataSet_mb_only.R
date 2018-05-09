#Make separate Cell Data Set for midbrain samples


##
## Merging the data
##
allCounts_mb <-  merge(mb_wt_1,mb_wt_2,by.x=0,by.y=0,all=TRUE)
allCounts_mb <-  merge(allCounts_mb,mb_ivs_1,by.x=1,by.y=0,all=TRUE)
allCounts_mb <-  merge(allCounts_mb,mb_ivs_2,by.x=1,by.y=0,all=TRUE)


##
## Writing the counts to file
##
rownames(allCounts_mb) <-  allCounts_mb[,1]
allCounts_mb <-  allCounts_mb[,-1]
write.csv(allCounts_mb,file='allCounts_mb.csv')

##
## Making the meta data
##
sData_mb <-  data.frame(
  'sampleName' = c('mb_wt_1','mb_wt_2','mb_ivs_1','mb_ivs_2'),
  'genotype' = c('wt','wt','ivs','ivs'),
  'region' = c(rep('mb',4)),
  'div' = rep(60,4),
  'genoReg' = c(rep('mbWt',2),(rep('mbIvs',2)))
)

allpData_mb <-  data.frame(
  'sampleName' = sub("-.*","",colnames(allCounts_mb)),
  'cellNum' = colnames(allCounts_mb)
)
allpData_mb <-  merge(allpData_mb,sData_mb,by.x=1,by.y=1,all.x=T,all.y=T,sort=F)
rownames(allpData_mb) = allpData_mb[,'cellNum']

## randomly sample the pData dataset
## to make sure that all the metadata looks good
## repeat this until you are convinced
allpData_mb[sample(1:dim(allpData_mb)[1],1),]


##
## filtering the data
library(Matrix)
## replace the NA values from the merge with 0
allCounts_mb[is.na(allCounts_mb)] = 0

## filtering out cells that do not have
## at least 100 reads.
validCells_mb <-  which(colSums(allCounts_mb) > 100)
length(validCells_mb)
##6396 cells total.
filteredCounts_mb <-  as(as.matrix(allCounts_mb[,names(validCells_mb)]), 'sparseMatrix')

## filtering out cells that do not fall within 2 sd of the
## mean number of reads
total_mRNAs_mb <-  Matrix::colSums(filteredCounts_mb)
upper_bound_mb <- 10^(mean(log10(total_mRNAs_mb))+
                        2*sd(log10(total_mRNAs_mb)))
lower_bound_mb <- 10^(mean(log10(total_mRNAs_mb)) -
                        2*sd(log10(total_mRNAs_mb)))

################ keeping only the cells in the upper and lower bound
filteredCounts_mb <-  filteredCounts_mb[,names(which(total_mRNAs_mb > lower_bound_mb & total_mRNAs_mb < upper_bound_mb))]

## filtering the genes
## requiring at least 4 read in at least
## 10 cells.
genesfData_mb <-  data.frame(rep("protein coding",dim(filteredCounts_mb)[1]),
                             row.names=rownames(filteredCounts_mb))
genesfData_mb[,'gene_short_name'] = rownames(filteredCounts_mb)

numCellExpressed_mb <-  Matrix::rowSums(filteredCounts_mb > 3)
genesfData_mb[,'numExpressed'] <-  numCellExpressed_mb
cellsPerGene_mb <-  Matrix::rowSums(filteredCounts_mb > 3)
expressedGenes_mb <-  names(which(cellsPerGene_mb > 10))
length(expressedGenes_mb)
#2801 genes expressed

##
## Loading into Monocle datframe
##
library(monocle)

##check that dimensions of the 3 input files are correct
dim(filteredCounts_mb[expressedGenes_mb,])
## 2801 rows 6144 columns
dim(allpData_mb[colnames(filteredCounts_mb),])
#should have tye same number of rows as Counts has.
#columns = 6144
dim(genesfData_mb[expressedGenes_mb,])
#should have the same number of rows as Counts has.
#rows =2801------------------continue running and editing here

HSMM_mb_1 <-  newCellDataSet(
  filteredCounts_mb[expressedGenes_mb,],
  phenoData = new("AnnotatedDataFrame",data=allpData_mb[colnames(filteredCounts_mb),]),
  featureData = new("AnnotatedDataFrame",data=genesfData_mb[expressedGenes_mb,]),
  expressionFamily=negbinomial())

HSMM_mb_1 <- detectGenes(HSMM_mb_1, min_expr = 4)
#sets global detection threshold to be used with this CellDataSet
HSMM_mb_1 <- estimateSizeFactors(HSMM_mb_1)
HSMM_mb_1 <- estimateDispersions(HSMM_mb_1)


qplot(total_mRNAs_mb[rownames(pData(HSMM_mb_1))], data=pData(HSMM_mb_1), color=genotype,geom="density")+
  geom_vline(xintercept=lower_bound_mb) +
  geom_vline(xintercept=upper_bound_mb)

qplot(total_mRNAs_mb[rownames(pData(HSMM_mb_1))], data=pData(HSMM_mb_1), color=sampleName,geom="density")+
  geom_vline(xintercept=lower_bound_mb) +
  geom_vline(xintercept=upper_bound_mb)




