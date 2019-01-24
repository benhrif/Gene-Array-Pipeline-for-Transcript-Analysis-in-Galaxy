
## load package

package.installed <- installed.packages()



suppressMessages(suppressWarnings(require(RSQLite)))
suppressMessages(suppressWarnings(require(DBI)))
suppressMessages(suppressWarnings(require(pvclust)))
suppressMessages(suppressWarnings(require(oligo)))
suppressMessages(suppressWarnings(require(Biobase)))


args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));

suppressMessages(suppressWarnings(load(args[1])));
Selected <- suppressWarnings(suppressMessages(as.matrix(read.table(args[2],header=TRUE))));
method_dist <- args[3]


if(exists("ppData")){
  dataQC <- ppData
  cluster.bootstrap <- suppressMessages(pvclust(exprs(dataQC), nboot=10, method.dist="abscor"))
}
if(exists("rawData")){
  dataQC <- rawData
  cluster.bootstrap <- suppressMessages(pvclust(exprs(dataQC), nboot=10, method.dist="abscor"))
}
if(exists("transcriptExpression")){
  dataQC <- transcriptExpression [Selected[,1],]
  
  cluster.bootstrap <- suppressMessages(pvclust(dataQC, nboot=10, method.dist=method_dist))
}

##output

pdf("Dendrogramm.pdf")

suppressMessages(suppressWarnings(plot(cluster.bootstrap, main="Cluster dendrogram of all samples")))
suppressMessages(suppressWarnings(pvrect(cluster.bootstrap)))

dev.off()
