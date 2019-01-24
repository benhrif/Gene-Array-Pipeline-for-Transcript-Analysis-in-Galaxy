#install packages
package.installed <- installed.packages()


#require library
suppressWarnings(suppressMessages(library("proxy")))
suppressWarnings(suppressMessages(library("methods")))
suppressWarnings(suppressMessages(library("S4Vectors")))
suppressWarnings(suppressMessages(library("oligo")))

# Get inputs
args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));
suppressWarnings(suppressMessages(load(args[1])));
Selected <- suppressWarnings(suppressMessages(as.matrix(read.table(args[2],header=TRUE))));
dist_methode <- args[3]

#heatmap

Expression_set <- exprs(ppData)

mat <- Expression_set[as.vector(Selected[,1]),]
main <- substr(args[2],1,nchar(args[2])-3)
main <- paste("Heatmap of",main,"using",dist_methode,"Distance")
pdf("heatmap.pdf")
heatmap(mat,dist(mat, method=dist_methode),main=paste("Heatmap using",dist_methode,"Distance"))
dev.off()