## Review of NCBI Listeria SNP clusters for David and Bill
## Errol Strain - 2015-10-01
## Read in the list of Listeria clusters

library(ape)
library(phangorn)

########################################################################################################
## Function to color clusters with food/env and clinical isolates that differ <= 20 SNPs
nodeColors <- function(cTree,distmat,metaDf) {
  unlist(sapply(cTree$tip.label,function(x) {
    ## Subset matrix to 
    temp <- distmat[distmat$target_acc_1==x | distmat$target_acc_2==x,]
    ## Subset 
    temp <- temp[temp$delta_positions_unambiguous<=20,]
    acc <- unique(c(temp$target_acc_1,temp$target_acc_2))
    
    if(length(grep("food",metaDf[metaDf$target_acc==x,"attribute_package"]))>0 | length(grep("Food",metaDf[metaDf$target_acc==x,"isolation_source"]))>0) {
      if(length(grep("clinical",metaDf[metaDf$target_acc %in% acc,"attribute_package"]))>0 | length(grep("Human",metaDf[metaDf$target_acc %in% acc,"isolation_source"]))>0) {
        return("blue")
      } else {
        return("black")
      }
    } else if (length(grep("clinical",metaDf[metaDf$target_acc==x,"attribute_package"]))>0 | length(grep("Human",metaDf[metaDf$target_acc==x,"isolation_source"]))>0) {
      if(length(grep("food",metaDf[metaDf$target_acc %in% acc,"attribute_package"]))>0 | length(grep("Food",metaDf[metaDf$target_acc %in% acc,"isolation_source"]))>0) {
        return("red")
      } else {
        return("black")
      }
    } else {
      return("black")
    }
  }))
}
########################################################################################################
metaFile <-list.files(pattern="metadata.tsv$")
bigMetaDf <- read.csv(metaFile[1],sep="\t",stringsAsFactors=FALSE)
dim(bigMetaDf)

distFile <-list.files(pattern="distances.tsv$")
distDf <- read.csv(distFile[1],sep="\t",stringsAsFactors=FALSE)
dim(distDf)

clusterList <- list.files(pattern="*.newick")
clusterList <- gsub(".newick_tree.newick","",clusterList)
numList <- as.numeric(gsub("PDS","",clusterList))
## Numerically sort the list
clusterList <- unlist(sapply(sort(numList), function(x) {return(clusterList[which(numList==x)])}))

pdf("NCBI_LmSNPclinfood.pdf")

for(i in 1:length(clusterList)) {
  treeName <- paste(clusterList[i],".newick_tree.newick",sep="")
  cTree <- read.tree(treeName)
  cTree$tip.label <- gsub("'","",unlist(cTree$tip.label))
  metaDf <- bigMetaDf[bigMetaDf$target_acc %in% cTree$tip.label,]

  if(cTree$Nnode > 1) { cTree <-midpoint(cTree) }
  distmat <- distDf[distDf$PDS_acc==clusterList[i],]
  cTree$tip.color <- nodeColors(cTree,distmat,metaDf)
  cTree$tip.label<- apply(metaDf[match(cTree$tip.label, metaDf$target_acc),c("target_acc","strain","collection_date","geo_loc_name","isolation_source")],1,paste,collapse = " ")
  plotTitle=clusterList[i]
  if("red" %in% cTree$tip.color && "blue" %in% cTree$tip.color) { 
    plotTitle=paste(plotTitle,"Clin-Food/Env")
    plot(cTree,main=plotTitle,tip.color=cTree$tip.color,cex=max(1-(nrow(metaDf)/30),0.2),align.tip.label=TRUE)
  }
}

dev.off()


