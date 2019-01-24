#install packages
package.installed <- installed.packages()

#require library
suppressWarnings(suppressMessages(library("VennDiagram")))

# Get inputs
args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));
Selected1 <- suppressWarnings(suppressMessages(as.matrix(read.table(args[1],header=TRUE))));
Selected2 <- suppressWarnings(suppressMessages(as.matrix(read.table(args[2],header=TRUE))));
l<-list(as.vector(Selected1[,1]),as.vector(Selected2[,1]))
names(l) <- c("List 1","List 2")
fill = c("yellow", "cornflowerblue")
cat_col = c("orange", "darkblue")
inter <- intersect(l[[1]],l[[2]])
if(args[3]!="None"){
	Selected3 <- suppressWarnings(suppressMessages(as.matrix(read.table(args[3],header=TRUE))));
	l<-list(l[[1]],l[[2]],as.vector(Selected3[,1]))
	names(l) <- c("List 1","List 2","List 3")
	inter <- intersect(l[[1]],l[[2]])
	inter <- intersect(inter,l[[3]])
fill = c("yellow", "cornflowerblue", "pink")
cat_col = c("orange", "darkblue", "pink4")
}
if(args[4]!="None"){
	Selected4 <- suppressWarnings(suppressMessages(as.matrix(read.table(args[4],header=TRUE))));
	l<-list(l[[1]],l[[2]],l[[3]],as.vector(Selected4[,1]))
	names(l) <- c("List 1","List 2","List 3","List 4")
	inter <- intersect(l[[1]],l[[2]])
	inter <- intersect(inter,l[[3]])
	inter <- intersect(inter,l[[4]])
fill = c("yellow", "cornflowerblue", "pink", "darkorchid1")
cat_col = c("orange", "darkblue", "pink4", "darkorchid4")
}
rownames(Selected1) <- Selected1[,1]
write.table(Selected1[inter,],"Intersection.txt", row.names=FALSE, sep="\t")
venn.diagram(l,"venn.png",col = "transparent",fill=fill,alpha = 0.5, cex = 1.5, fontfamily = "serif", 
             cat.fontface = "bold", cat.col = cat_col, cat.cex = 1.5, cat.pos = 0, cat.dist = 0.03,
             cat.fontfamily = "serif", imagetype="png", main="Venn Diagram", main.pos=c(0.5,1.05),
             rotation.degree = 270,
             margin = 0.2
			 )