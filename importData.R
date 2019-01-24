#folder <- "/home/anygenes/Desktop/affymetrix"
## load package
package.installed <- installed.packages()

suppressWarnings(suppressMessages(library("oligo")))
suppressWarnings(suppressMessages(library("S4Vectors")))
args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));

folder <- args[1];

Files<-suppressMessages(suppressWarnings(list.celfiles(folder,full.name=TRUE)))

rawData<-suppressMessages(suppressWarnings(read.celfiles(Files)))
colnames(rawData) <- substr(colnames(rawData),1,nchar(colnames(rawData))-4)
suppressMessages(suppressWarnings(save(rawData,file="rawData.rdata")))
