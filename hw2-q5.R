
# clear the environment
rm(list = ls())

# install the packages
if (F) {
  install.packages("GEOquery")
}

#---------------------#
# Question 1
#
#---------------------#

library(GEOquery)

# load the data
# Berry UK Training (array)
gse19439 = getGEO("GSE19439",GSEMatrix = T)
show(gse19439)
pdat19439 = pData(phenoData(gse19439[[1]]))[,c(1,11:16)]
#---------Notation-------------#
# 1=title, 11=gender(Female, Male), 12=ethinicity(White,Black,South Asian,Asian Other,Other),
# 13=illness(Control (BCG-),Control (BCG+), Latent,PTB), 14=geographical region(London),
# 15=bcg vaccinated(No,Yes,Not_Known), 16=region of birth(UK_Born,non_UK_Born)
#------------------------------#
colnames(pdat19439) <- c("group", "gender", "ethnicity", "illness", "geographical_region", "bcg_vaccinated", "region_of_birth")








#------------------------#
# Question 5
#------------------------#
rm(list = ls())
setwd("C:/Users/amanda/Dropbox/Group2/HW2/JZhang")
load("HW2X.RData")
colnames(X) = as.character(c(1:dim(X)[2])) # name the columns
rownames(X)=as.character(c(1:dim(X)[1])) # name the rows
heatmap(X,Rowv=NA,Colv=NA,hclustfun = hclust) # total mess
heatmap(X,Colv=NA,hclustfun = hclust) # reorder the rows only # still mess
heatmap(X,Rowv=NA,hclustfun = hclust) # reorder the columns only

h1=hclust(dist(X))
plot(h1)
groups = cutree(h1, k=2) # cut tree into 4 clusters
which(groups==1) # 1:30 & 211:240



# the bottom and top 30 row belong to the same clust
# So let's trying to make all those block smooth.


# order the matrix by col means of first 30 rows
colmean = colMeans(X[c(1:30),])
XMean = X[,order(colmean)]
heatmap(XMean,Rowv = NA,Colv = NA)

# use correlation instead of means to catch more information.
# Choose the first column as beginning
# which is the column with minimun mean

colnames(XMean)[1] # 73 column in X
n = dim(X)[2] # the number of columns
orderCol_1 =rep(0,n)
orderCol_1[1] = 73
remainingCol = X[,-73]
currentCol = X[,73]
for (i in 1:(n - 1))
{
  corCol = 0
  #get all corrlation
  for (j in 1:(n - i))
  {
    corCol[j] = cor(currentCol,remainingCol[,j],method = "pearson")
  }
  maxInd = which.max(corCol)
  
  currentCol = remainingCol[,maxInd,drop = F]
  remainingCol = remainingCol[,-maxInd,drop = F]
  orderCol_1[i + 1] = as.numeric(colnames(currentCol))
  
}
heatmap(X[,orderCol_1],Rowv = NA,Colv = NA)




# Also, we can choose the last column as the beginning
# which is the column with maximum mean
# should be symmetric
colnames(XMean)[140] # 9 column in X
n = dim(X)[2] # the number of columns
orderCol_2 =rep(0,n)
orderCol_2[1] = 9
remainingCol = X[,-9]
currentCol = X[,9]
for (i in 1:(n - 1))
{
  corCol = 0
  #get all corrlation
  for (j in 1:(n - i))
  {
    corCol[j] = cor(currentCol,remainingCol[,j],method = "pearson")
  }
  maxInd = which.max(corCol)
  
  currentCol = remainingCol[,maxInd,drop = F]
  remainingCol = remainingCol[,-maxInd,drop = F]
  orderCol_2[i + 1] = as.numeric(colnames(currentCol))
  
}
heatmap(X[,orderCol_2],Rowv = NA,Colv = NA)


