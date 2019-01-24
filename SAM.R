# Get inputs
package.installed <- installed.packages()

suppressWarnings(suppressMessages(library("siggenes")))
suppressWarnings(suppressMessages(library("methods")))
suppressWarnings(suppressMessages(library("S4Vectors")))
args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));
suppressWarnings(suppressMessages(load(args[1])));
Design_matrix <- suppressWarnings(suppressMessages(as.matrix(read.table(args[2],header=TRUE))));
delta <- suppressWarnings(suppressMessages(as.numeric(args[3])))
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

#require library
suppressWarnings(suppressMessages(source(paste(args[4],"/annoter_gene.R",sep=""))))
suppressWarnings(suppressMessages(source(paste(args[4],"/annoter_variants.R",sep=""))))
#SAM Function Arguments
if(exists("ppData")){
  Expression_set <- exprs(ppData)
}
if(exists("transcriptExpression")){
  Expression_set <- transcriptExpression
}
rownames(Design_matrix) <- Design_matrix[,1]
col <- colnames(Expression_set)
cl <- vector()
listControl <- Design_matrix[which(Design_matrix[col,2]=="C"),1]
listTrait <- Design_matrix[which(Design_matrix[col,2]=="T"),1]
cl <- c(rep(0,length(listControl)),rep(1,length(listTrait)))

#SAM Matrix
data_sam <-Expression_set[,c(listControl,listTrait)]

#Test SAM
sam_output <- suppressWarnings(suppressMessages(sam(data_sam, cl, control = samControl (n.delta = 6))))

sam_result <- cbind(P.value=sam_output@p.value, FC=sam_output@fold)
s1 <-summary(sam_output,delta=delta)
listeSam <- names(s1@row.sig.genes)

SelectedList <- sam_result[listeSam,]

if(exists("ppData")){
  listeFinalSAM <- annoter_gene(SelectedList,Database)
}

#FC
result <- Expression_set[listeFinalSAM[,1],]

FC_result <- vector()
for(i in 1:dim(result)[1]){
  if(mean(result[i,listControl]) >= mean(result[i,listTrait]))
    FC_result[i] <- 2^(mean(result[i,listControl])-mean(result[i,listTrait]))
  else
    FC_result[i] <- -1/(2^(mean(result[i,listControl])-mean(result[i,listTrait])))
}
listeFinalSAM[,3] <- FC_result
listeFinalSAM <- listeFinalSAM[order(abs(listeFinalSAM[,3]),decreasing = TRUE),]
colnames(listeFinalSAM)<-c("Gene ID" ,"P.Value", "Fold Change", "Symbol","ACCNUM","DESC","ENSEMBL")
write.table(listeFinalSAM,"SelectedListSAM.txt",row.names = FALSE, sep="\t")

pdf("SAMPlot.pdf")

sam.plot2(sam_output,delta,main=paste("SAM Plot for delta =",delta,sep=" "),xlab="Expected values ", ylab="Observed Values ")

dev.off()
