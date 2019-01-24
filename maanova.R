#install packages
package.installed <- installed.packages()

# Get inputs

suppressWarnings(suppressMessages(library("qvalue")))
suppressWarnings(suppressMessages(library("maanova")))
suppressWarnings(suppressMessages(library("methods")))
suppressWarnings(suppressMessages(library("S4Vectors")))
args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));
suppressWarnings(suppressMessages(load(args[1])));
Design_matrix <- suppressWarnings(suppressMessages(as.matrix(read.table(args[2],header=TRUE))));
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

#require functions
suppressWarnings(suppressMessages(source(paste(args[5],"/annoter_gene.R",sep=""))))
suppressWarnings(suppressMessages(source(paste(args[5],"/annoter_variants.R",sep=""))))

Pvalue <- as.numeric(args[3])
FoldChange <- as.numeric(args[4])

rownames(Design_matrix) <- Design_matrix[,1]
if(exists("ppData")){
Expression_set <- exprs(ppData)
}
col <- colnames(Expression_set)
cl <- vector()
listControl <- as.vector(Design_matrix[which(Design_matrix[col,2]=="C"),1])
listTrait <- as.vector(Design_matrix[which(Design_matrix[col,2]=="T"),1])

cl <- c(rep(0,length(listControl)),rep(1,length(listTrait)))
Expression_set <- Expression_set[,c(listControl,listTrait)]
##########################################################################################################################################
#MAANOVA
##########################################################################################################################################
#design matrix
write.table(Expression_set,paste(args[5],"/normalized_data.dat",sep=""),sep="\t",row=F,quote=F)
groups <- Design_matrix[,2]
design.matrix <- data.frame(Array=Design_matrix[,1], Strain=Design_matrix[,2])
write.table(design.matrix,paste(args[5],"/design.dat",sep=""),sep="\t",row=F,quote=F)
norm_data <- read.madata(paste(args[5],"/normalized_data.dat",sep=""),designfile=paste(args[5],"/design.dat",sep=""),probeid=row.names(Expression_set),
                         intensity=1,cloneid=1,pmt=2,spot=F)
#mannova
fit.fix <- suppressWarnings(suppressMessages(fitmaanova(norm_data, formula = ~Strain)))
ftest.all <- suppressWarnings(suppressMessages(matest(norm_data, fit.fix, term="Strain", n.perm= 50)))
test.fix <- adjPval(ftest.all, method = "jsFDR")
#summary result
summary_<-suppressWarnings(suppressMessages(summarytable(ftest.all,outfile=paste(args[5],"/summarytable.csv",sep=""))))
#volcano plot
pdf("VolcanotPlotManova.pdf")

idx <- volcano(test.fix, title="volcano plot",threshold=c(Pvalue,Pvalue),method = "unadj")
abline(v=c(-FoldChange,FoldChange),col="red")

dev.off()

#select genes from volcano where p-val > Pvalue
pval_mannova <- getPval.volcano(test.fix, "unadj", 1)
rownames(pval_mannova) <- rownames(Expression_set)
IDX <- idx$idx.all
resu <- row.names(Expression_set[IDX,])
st <- cbind("p-value"=pval_mannova[resu,],"Fold change"=summary_[resu,1])
rownames(st)<-resu
# selected genes where abs(fc) > Fold Change
st_selescted_fc<-st[abs(st[,2])>=FoldChange,]

if(exists("ppData")){
  listFinalManova <- annoter_gene(st_selescted_fc,Database)
  listFinalManova <- listFinalManova[order(abs(listFinalManova[,3]),decreasing = TRUE),]
	colnames(listFinalManova)<-c("Gene ID" ,"P.Value", "Fold Change", "Symbol","ACCNUM","DESC","ENSEMBL")
}
write.table(listFinalManova,"SelectedListManova.txt",row.names = FALSE, sep="\t")

