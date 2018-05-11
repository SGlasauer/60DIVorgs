###classifying cells with known markers

#set up expression ids
PAX6_id <- row.names(subset(fData(HSMM1), gene_short_name == "PAX6"))
NES_id <- row.names(subset(fData(HSMM1), gene_short_name == "NES"))


#register a set of classification functions
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "RadialGlia1",
                   classify_func = function(x) {
                     x[PAX6_id,] >= 1})

cth <- addCellType(cth, "RadialGlia2",
                   classify_func = function(x) {
                     x[NES_id,] >=1})

#Classify each cell in the CellDataSet according to these types
HSMM1 <- classifyCells(HSMM1, cth, 0.1)

#check how many cells are in each class
table(pData(HSMM1)$CellType)
#ambiguous 247, RadialGlia1 590, RadialGlia2 2128, Unknown 9283

#plot as a pie chart
pie <- ggplot(pData(HSMM1),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


#######################################################################
#######################################################################
##now do the same for forebrain

PAX6_id_fb <- row.names(subset(fData(HSMM_fb_1), gene_short_name == "PAX6"))
NES_id_fb <- row.names(subset(fData(HSMM_fb_1), gene_short_name == "NES"))


#register a set of classification functions
cth_fb <- newCellTypeHierarchy()
cth_fb <- addCellType(cth_fb, "RadialGlia1_fb",
                   classify_func = function(x) {
                     x[PAX6_id_fb,] >= 1})

cth_fb <- addCellType(cth_fb, "RadialGlia2_fb",
                   classify_func = function(x) {
                     x[NES_id_fb,] >=1})

#Classify each cell in the CellDataSet according to these types
HSMM_fb_1 <- classifyCells(HSMM_fb_1, cth_fb, 0.1)

#check how many cells are in each class
table(pData(HSMM_fb_1)$CellType)
#ambiguous 215, RadialGlia1 552, RadialGlia2 1187, Unknown 4149

#plot as a pie chart
pie_fb <- ggplot(pData(HSMM_fb_1),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie_fb + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


#######################################################################
#######################################################################
##now do the same for midbrain

PAX6_id_mb <- row.names(subset(fData(HSMM_mb_1), gene_short_name == "PAX6"))
NES_id_mb <- row.names(subset(fData(HSMM_mb_1), gene_short_name == "NES"))


#register a set of classification functions
cth_mb <- newCellTypeHierarchy()
cth_mb <- addCellType(cth_mb, "RadialGlia1_mb",
                      classify_func = function(x) {
                        x[PAX6_id_mb,] >= 1})

cth_mb <- addCellType(cth_mb, "RadialGlia2_mb",
                      classify_func = function(x) {
                        x[NES_id_mb,] >=1})

#Classify each cell in the CellDataSet according to these types---------error message!
HSMM_mb_1 <- classifyCells(HSMM_mb_1, cth_mb, 0.1)

#check how many cells are in each class
table(pData(HSMM_mb_1)$CellType)


#plot as a pie chart
pie_mb <- ggplot(pData(HSMM_mb_1),
                 aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie_mb + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())




#group the cells by the pData table column "Cluster". Apply classification funcitons to the cell clusters.
#If a group is at least 5% of a type, make them all that type.
#If the group is 5% one type and 5% one other, mark the cluster "Ambiguous"
#DO NOT RUN NOW - only if clusters are assigned
HSMM1 <- classifyCells(mix, Cluster, 0.05)
