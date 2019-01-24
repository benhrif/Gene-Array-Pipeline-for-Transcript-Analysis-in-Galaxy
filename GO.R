#install packages
package.installed <- installed.packages()



suppressWarnings(suppressMessages(library("topGO")))
suppressWarnings(suppressMessages(library("oligo")))
suppressWarnings(suppressMessages(library("Rgraphviz")))
suppressWarnings(suppressMessages(library("lattice")))

args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));

suppressWarnings(suppressMessages(load(args[1])));
Selected <- suppressWarnings(suppressMessages(as.matrix(read.table(args[2],header=TRUE))));
Design_matrix <- suppressWarnings(suppressMessages(as.matrix(read.table(args[3],header=TRUE))));
suppressWarnings(suppressMessages(source(paste(args[4],"/barchart.R",sep=""))))
suppressWarnings(suppressMessages(source(paste(args[4],"/printGraphGO.R",sep=""))))
eset <- exprs(ppData)
eset <- eset[as.vector(Selected[,1]),]
if(exists("ppData")){
  chip <- gsub('\\.','',ppData@annotation)
  chip <- substr(chip,3,nchar(chip))
  Database <- paste(chip,"transcriptcluster.db",sep="")
  if(!Database %in% rownames(package.installed)){
    suppressMessages(suppressWarnings(source("http://bioconductor.org/biocLite.R")))
    suppressMessages(suppressWarnings(biocLite(Database)))
  }
  suppressMessages(library(Database,character.only = TRUE))
}

affyLib <- Database
rownames(Design_matrix) <- Design_matrix[,1]
col <- colnames(eset)
y <- vector()
listControl <- Design_matrix[which(Design_matrix[col,2]=="C"),1]
listTrait <- Design_matrix[which(Design_matrix[col,2]=="T"),1]
y <- c(rep(0,length(listControl)),rep(1,length(listTrait)))

geneList <-     suppressMessages(suppressWarnings(getPvalues(eset, classlabel = y,test="f")))

graphBP <- suppressMessages(suppressWarnings(printGraphGO(geneList,"BP",Database,"BP")))
graphMF <- suppressMessages(suppressWarnings(printGraphGO(geneList,"MF",Database,"MF")))
graphCC <- suppressMessages(suppressWarnings(printGraphGO(geneList,"CC",Database,"CC")))


pdf("Barchart_BiologicalProcesses.pdf")

gt_genesBP <- suppressMessages(suppressWarnings(genes_barchart(graphBP[[1]],graphBP[[2]],"BP",Database)))

dev.off()

pdf("Barchart_MoleculaireFunction.pdf")

gt_genesMF <- suppressMessages(suppressWarnings(genes_barchart(graphMF[[1]],graphMF[[2]],"MF",Database)))

dev.off()

pdf("Barchart_CellulaireComponent.pdf")

gt_genesCC <- suppressMessages(suppressWarnings(genes_barchart(graphCC[[1]],graphCC[[2]],"CC",Database)))

dev.off()
