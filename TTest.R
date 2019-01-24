# Get inputs
package.installed <- installed.packages()

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
suppressWarnings(suppressMessages(source(paste(args[5],"/annoter_gene.R",sep=""))))
suppressWarnings(suppressMessages(source(paste(args[5],"/annoter_variants.R",sep=""))))

Pvalue <- as.numeric(args[3])
FoldChange <- as.numeric(args[4])

rownames(Design_matrix) <- Design_matrix[,1]

Expression_set <- exprs(ppData)

col <- colnames(Expression_set)
cl <- vector()
listControl <- Design_matrix[which(Design_matrix[col,2]=="C"),1]
listTrait <- Design_matrix[which(Design_matrix[col,2]=="T"),1]

cl <- c(rep(0,length(listControl)),rep(1,length(listTrait)))

n<-dim(Expression_set)[1]
FC <- matrix(c(rep(0, n)), ncol = 1)
p.value <- matrix(c(rep(0, n)), ncol = 1)
stat <- matrix(c(rep(0, n)), ncol = 1)
for (i in 1:n) {
  x1 <- Expression_set[i,listControl]
  x2 <- Expression_set[i, listTrait]
  if ((sd(x1) + sd(x2)) == 0) {
    p.value[i] <- 1
    stat[i] <- 0
  }else {
    tt <- t.test(x1, x2, var.equal = TRUE)
    p.value[i] <- tt$p.value
    stat[i] <- tt$statistic
  }
  FC[i] <- 2^(mean(x2))/2^(mean(x1))
}

FTT<-list(idnames = rownames(Expression_set), FC = log2(FC), p.value = p.value, stat = stat)

FTT<-data.frame(FTT)

liste_gene <-FTT[FTT$p.value<Pvalue & abs(FTT$FC)>log2(FoldChange),]
rownames(liste_gene) <- liste_gene[,1]
SelectedList <- liste_gene[,-1]

  listeFinalTtest <- annoter_gene(SelectedList,Database)

listeFinalTtest <- listeFinalTtest[which(listeFinalTtest[,5]!="NA"),]
listeFinalTtest <- listeFinalTtest[,-4]
listeFinalTtest <- listeFinalTtest[order(abs(listeFinalTtest[,2]),decreasing = TRUE),]
colnames(listeFinalTtest)<-c("Gene ID" ,"Fold Change","P.Value", "Symbol","ACCNUM","DESC","ENSEMBL")
write.table(listeFinalTtest,"SelectedListTtest.txt",row.names = FALSE, sep="\t")

pdf("VolcanotPlot.pdf")

volcano1 <- plot(FTT$FC, -log10(FTT$p.value), main="Ttest Volcano Plot" ,xlab="log2 fold change", ylab="-log10 p-value", col=(ifelse(abs(FTT$FC)>log2(FoldChange) & FTT$p.value<Pvalue ,"green","black")))
abline(v=c(-log2(FoldChange),log2(FoldChange)), col="purple")
abline(h=-log10(Pvalue), col="red")

dev.off()