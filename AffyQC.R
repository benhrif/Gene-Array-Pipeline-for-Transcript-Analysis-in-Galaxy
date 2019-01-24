
## load package
suppressMessages(suppressWarnings(require(RSQLite)))
suppressMessages(suppressWarnings(require(DBI)))
suppressMessages(suppressWarnings(require(pvclust)))
suppressMessages(suppressWarnings(require(oligo)))
suppressMessages(suppressWarnings(require(Biobase)))


args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));

suppressMessages(suppressWarnings(load(args[1])));
suppressMessages(suppressWarnings(load(args[2])));

if(exists("rawData")){
   dataQC1 <- rawData
}
if(exists("ppData")){
   dataQC2 <- ppData
}
##output

pdf("AffyQCPlots.pdf")
suppressMessages(suppressWarnings(boxplot(dataQC1, main="Box Plot of all samples Before Preprocessing",las=2)))
suppressMessages(suppressWarnings(boxplot(dataQC2, main="Box Plot of all samples After Preprocessing",las=2)))
dev.off()
