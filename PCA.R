#install packages
package.installed <- installed.packages()



#require library
suppressWarnings(suppressMessages(library("proxy")))
suppressWarnings(suppressMessages(library("FactoMineR")))
suppressWarnings(suppressMessages(library("gplots")))
suppressWarnings(suppressMessages(library("oligo")))

# Get inputs
args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));
suppressWarnings(suppressMessages(load(args[1])));
Selected <- suppressWarnings(suppressMessages(as.matrix(read.table(args[2],header=TRUE))));

#PCA
Expression_set <- exprs(ppData)
mat <- Expression_set[as.vector(Selected[,1]),]

pca_liste <- suppressWarnings(suppressMessages(PCA(mat,scale.unit=TRUE, ncp=5, graph=T)))

pdf("pca.pdf")
plot.PCA(pca_liste, choix="var", cex=1, new.plot=FALSE)
plot.PCA(pca_liste, choix="ind", cex=0.5, new.plot=FALSE)
barplot(pca_liste$eig[,2], names=paste("Dim",1:nrow(pca_liste$eig)))
dev.off()
