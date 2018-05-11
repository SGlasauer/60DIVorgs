#This script is for reading in expression files, and preparing the CellDataSet for analysis in Monocle

##pwd
##"/data/Lisi/RNAseqData/StellasOrganoids/60DIV/Analysis"

##/bin/R

##
## Reading in the files
##
fb_wt_1 <-  read.table("/data/Lisi/RNAseqData/StellasOrganoids/60DIV/60DIVorganoids/FB_WT_1/star_gene_exon_tagged.dge.txt.gz",header=T)
rownames(fb_wt_1) =  fb_wt_1[,1]
fb_wt_1 =  fb_wt_1[,-1]
colnames(fb_wt_1) =  paste("fb_wt_1-",1:(dim(fb_wt_1)[2]),sep="")

fb_wt_2 <-  read.table("/data/Lisi/RNAseqData/StellasOrganoids/60DIV/60DIVorganoids/FB_WT_2/star_gene_exon_tagged.dge.txt.gz",header=T)
rownames(fb_wt_2) =  fb_wt_2[,1]
fb_wt_2 =  fb_wt_2[,-1]
colnames(fb_wt_2) =  paste("fb_wt_2-",1:(dim(fb_wt_2)[2]),sep="")

fb_ivs_1 <-  read.table("/data/Lisi/RNAseqData/StellasOrganoids/60DIV/60DIVorganoids/FB_IVS_1/star_gene_exon_tagged.dge.txt.gz",header=T)
rownames(fb_ivs_1) =  fb_ivs_1[,1]
fb_ivs_1 =  fb_ivs_1[,-1]
colnames(fb_ivs_1) =  paste("fb_ivs_1-",1:(dim(fb_ivs_1)[2]),sep="")

fb_ivs_2 <-  read.table("/data/Lisi/RNAseqData/StellasOrganoids/60DIV/60DIVorganoids/FB_IVS_2/star_gene_exon_tagged.dge.txt.gz",header=T)
rownames(fb_ivs_2) =  fb_ivs_2[,1]
fb_ivs_2 =  fb_ivs_2[,-1]
colnames(fb_ivs_2) =  paste("fb_ivs_2-",1:(dim(fb_ivs_2)[2]),sep="")

mb_wt_1 <-  read.table("/data/Lisi/RNAseqData/StellasOrganoids/60DIV/60DIVorganoids/MB_WT_1/star_gene_exon_tagged.dge.txt.gz",header=T)
rownames(mb_wt_1) =  mb_wt_1[,1]
mb_wt_1 =  mb_wt_1[,-1]
colnames(mb_wt_1) =  paste("mb_wt_1-",1:(dim(mb_wt_1)[2]),sep="")

mb_wt_2 <-  read.table("/data/Lisi/RNAseqData/StellasOrganoids/60DIV/60DIVorganoids/MB_WT_2/star_gene_exon_tagged.dge.txt.gz",header=T)
rownames(mb_wt_2) <-  mb_wt_2[,1]
mb_wt_2 = mb_wt_2[,-1]
colnames(mb_wt_2) <-  paste("mb_wt_2-",1:(dim(mb_wt_2)[2]),sep="")

mb_ivs_1 <-  read.table("/data/Lisi/RNAseqData/StellasOrganoids/60DIV/60DIVorganoids/MB_IVS_1/star_gene_exon_tagged.dge.txt.gz",header=T)
rownames(mb_ivs_1) =  mb_ivs_1[,1]
mb_ivs_1 = mb_ivs_1[,-1]
colnames(mb_ivs_1) =  paste("mb_ivs_1-",1:(dim(mb_ivs_1)[2]),sep="")

mb_ivs_2 <-  read.table("/data/Lisi/RNAseqData/StellasOrganoids/60DIV/60DIVorganoids/MB_IVS_2/star_gene_exon_tagged.dge.txt.gz",header=T)
rownames(mb_ivs_2) =  mb_ivs_2[,1]
mb_ivs_2 = mb_ivs_2[,-1]
colnames(mb_ivs_2) =  paste("mb_ivs_2-",1:(dim(mb_ivs_2)[2]),sep="")


##
## Merging the data
##
allCounts <-  merge(fb_wt_1,fb_wt_2,by.x=0,by.y=0,all=TRUE)
allCounts <-  merge(allCounts,fb_ivs_1,by.x=1,by.y=0,all=TRUE)
allCounts <-  merge(allCounts,fb_ivs_2,by.x=1,by.y=0,all=TRUE)
allCounts <-  merge(allCounts,mb_wt_1,by.x=1,by.y=0,all=TRUE)
allCounts <-  merge(allCounts,mb_wt_2,by.x=1,by.y=0,all=TRUE)
allCounts <-  merge(allCounts,mb_ivs_1,by.x=1,by.y=0,all=TRUE)
allCounts <-  merge(allCounts,mb_ivs_2,by.x=1,by.y=0,all=TRUE)

##
## Writing the counts to file
##
rownames(allCounts) <-  allCounts[,1]
allCounts <-  allCounts[,-1]
write.csv(allCounts,file='allCounts.csv')

##
## Making the meta data
##
sData <-  data.frame(
  'sampleName' = c('fb_wt_1','fb_wt_2','fb_ivs_1','fb_ivs_2',
                   'mb_wt_1','mb_wt_2','mb_ivs_1','mb_ivs_2'),
  'genotype' = rep(c('wt','wt','ivs','ivs'),2),
  'region' = c(rep('fb',4),rep('mb',4)),
  'div' = rep(60,8),
  'genoReg' = c(rep('fbWt',2), rep ('fbIvs', 2), rep ('mbWt', 2), rep ('mbIVS', 2))
  )

allpData <-  data.frame(
  'sampleName' = sub("-.*","",colnames(allCounts)),
  'cellNum' = colnames(allCounts)
)
allpData <-  merge(allpData,sData,by.x=1,by.y=1,all.x=T,all.y=T,sort=F)
rownames(allpData) = allpData[,'cellNum']

## randomly sample the pData dataset
## to make sure that all the metadata looks good
## repeat this until you are convinced
allpData[sample(1:dim(allpData)[1],1),]


##
## filtering the data
library(Matrix)
## replace the NA values from the merge with 0
allCounts[is.na(allCounts)] = 0

## filtering out cells that do not have
## at least 100 reads.
validCells <-  which(Matrix::colSums(allCounts) > 100)
length(validCells)
##12793 cells total.
filteredCounts <-  as(as.matrix(allCounts[,names(validCells)]),'sparseMatrix')

## filtering out cells that do not fall within 2 sd of the
## mean number of reads
total_mRNAs <-  Matrix::colSums(filteredCounts)
upper_bound <- 10^(mean(log10(total_mRNAs))+
                     2*sd(log10(total_mRNAs)))
lower_bound <- 10^(mean(log10(total_mRNAs)) -
                     2*sd(log10(total_mRNAs)))

################ keeping only the cells in the upper and lower bound
filteredCounts = filteredCounts[,names(which(total_mRNAs > lower_bound & total_mRNAs < upper_bound))]

## filtering the genes
## requiring at least 4 read in at least
## 10 cells.
genesfData <-  data.frame(rep("protein coding",dim(filteredCounts)[1]),
                        row.names=rownames(filteredCounts))
genesfData[,'gene_short_name'] = rownames(filteredCounts)

numCellExpressed <-  Matrix::rowSums(filteredCounts > 3)
genesfData[,'numExpressed'] <-  numCellExpressed
cellsPerGene <-  Matrix::rowSums(filteredCounts > 3)
expressedGenes <-  names(which(cellsPerGene > 10))
length(expressedGenes)
#3648 genes expressed

##
## Loading into Monocle datframe
##
library(monocle)

##check that dimensions of the 3 input files are correct
dim(filteredCounts[expressedGenes,])
## 3648 rows 12248 columns
dim(allpData[colnames(filteredCounts),])
#should have tye same number of rows as Counts has.
#columns = 12248
dim(genesfData[expressedGenes,])
#should have the same number of rows as Counts has.
#rows =3648

HSMM1 <-  newCellDataSet(
  filteredCounts[expressedGenes,],
  phenoData = new("AnnotatedDataFrame",data=allpData[colnames(filteredCounts),]),
  featureData = new("AnnotatedDataFrame",data=genesfData[expressedGenes,]),
  expressionFamily=negbinomial())

HSMM1 <- detectGenes(HSMM1, min_expr = 4)
HSMM1 <- estimateSizeFactors(HSMM1)
HSMM1 <- estimateDispersions(HSMM1)

qplot(total_mRNAs[rownames(pData(HSMM1))], data=pData(HSMM1), color=genotype,geom="density")+
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound)

qplot(total_mRNAs[rownames(pData(HSMM1))], data=pData(HSMM1), color=genoReg,geom="density")+
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound)

qplot(total_mRNAs[rownames(pData(HSMM1))], data=pData(HSMM1), color=sampleName,geom="density")+
  geom_vline(xintercept=lower_bound) +
  geom_vline(xintercept=upper_bound)

