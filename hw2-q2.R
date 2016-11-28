## QUESTION 2
load("gse19439.RData")
load("pdat19439.RData")
assay.data <- exprs(gse19439[[1]])

feature.selection <- function(feature, subgroup){assay.data[,colnames(assay.data) %in% rownames(pdat19439[feature == subgroup,])]}

cluster.analysis <- function(featurelabel, distance.metric, cluster.metric){
        colnames(assay.data) <- paste(colnames(assay.data), pdat19439[[featurelabel]], sep="-")
        dissimilarity <- 1-cor(assay.data)
        dist <- dist(dissimilarity, method = distance.metric)
        plot(hclust(dist, method=cluster.metric), main=NULL)
        title(paste0(distance.metric, " distance dendrogram (", cluster.metric, " linkage) labeled by ", featurelabel))
}

cluster.analysis("illness", "manhattan", "complete")

cluster.analysis.subset <- function(feature, feature.subgroup, featurelabel, distance.metric, cluster.metric){
        select <- feature.selection(feature, feature.subgroup)
        new.pdat <- pdat19439[rownames(pdat19439) %in% colnames(select),]
        colnames(select) <- paste(colnames(select), new.pdat[[featurelabel]], sep="-")
        dissimilarity <- 1-cor(select)
        dist <- dist(dissimilarity, method = distance.metric)
        plot(hclust(dist, method=cluster.metric), main=NULL)
        title(paste0(distance.metric, " distance dendrogram (", cluster.metric, " linkage) subsetted by ", feature.subgroup, " and labeled by ", featurelabel))
}

cluster.analysis.subset(pdat19439$illness, "PTB", "ethnicity", "euclidean", "complete")

