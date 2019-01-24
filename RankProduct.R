#install packages
package.installed <- installed.packages()

# Get inputs
suppressWarnings(suppressMessages(library("methods")))
suppressWarnings(suppressMessages(library("S4Vectors")))
suppressWarnings(suppressMessages(library("RankProd")))
args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));
suppressWarnings(suppressMessages(load(args[1])));
Design_matrix <- suppressWarnings(suppressMessages(as.matrix(read.table(args[2],header=TRUE))));
PFP <- suppressWarnings(suppressMessages(as.numeric(args[3])))
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

#Rank Product Function Arguments
if(exists("ppData")){
  Expression_set <- exprs(ppData)
}

rownames(Design_matrix) <- Design_matrix[,1]
col <- colnames(Expression_set)
cl <- vector()
listControl <- as.vector(Design_matrix[which(Design_matrix[col,2]=="C"),1])
listTrait <- as.vector(Design_matrix[which(Design_matrix[col,2]=="T"),1])
cl <- c(rep(0,length(listControl)),rep(1,length(listTrait)))
Expression_set <- Expression_set[,c(listControl,listTrait)]
##########################################################################################################################################
#Rank Product
##########################################################################################################################################
##RP
rankProd_output  <- RP (Expression_set, cl, num.perm=10, rand=123, logged = F,plot=FALSE, gene.names = rownames(Expression_set))
##top genes with PFP
tables <- topGene(rankProd_output,cutoff=PFP,method="pfp", logged= F, gene.names = rownames(Expression_set))
#results files
tables$Table1[,1]<- rownames(tables$Table1)
colnames(tables$Table1)<-c("gene ID" ,"RP/Rsum", "FC", "pfp", "p.value")
tables$Table2[,1]<- rownames(tables$Table2)
colnames(tables$Table2)<-c("gene ID" ,"RP/Rsum", "FC", "pfp", "p.value")
#up- and down-file

##plot
pdf(file="RankProductPlot.pdf")
plotRP(rankProd_output, cutoff=PFP)
dev.off()
#up regulated + down regulated
SelectedList <- rbind(tables$Table1,tables$Table2)
if(exists("ppData")){
  listeFinalRP <- annoter_gene(SelectedList,Database)
  listeFinalRP <- listeFinalRP[,-c(1,3,5)]
  colnames(listeFinalRP)<-c("Gene ID" ,"Fold Change", "pfp", "Symbol","ACCNUM","DESC","ENSEMBL")
}
###########
#FC
result <- Expression_set[listeFinalRP[,1],]

FC_result <- vector()
for(i in 1:dim(result)[1]){
  if(mean(result[i,listControl]) >= mean(result[i,listTrait]))
    FC_result[i] <- 2^(mean(result[i,listControl])-mean(result[i,listTrait]))
  else
    FC_result[i] <- -1/(2^(mean(result[i,listControl])-mean(result[i,listTrait])))
}
listeFinalRP[,2] <- FC_result
listeFinalRP <- listeFinalRP[order(abs(listeFinalRP[,2]),decreasing = TRUE),]
colnames(listeFinalRP)<-c("Gene ID" ,"Fold Change","PFP", "Symbol","ACCNUM","DESC","ENSEMBL")	  
write.table(listeFinalRP,"SelectedListRankProduct.txt", row.names = FALSE, sep="\t")





