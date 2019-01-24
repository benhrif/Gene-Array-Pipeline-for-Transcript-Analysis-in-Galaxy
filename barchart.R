genes_barchart<-function(allRes,sampleGOdata,GO,Database){
  
  affyLib<-Database
  gt_total <- data.frame()
  gt_genes <- list()
  for(i in 1:dim(allRes)[1]){
    gt <- printGenes(sampleGOdata,
                       whichTerms = allRes[i, "GO.ID"], 
                       chip = affyLib, numChar = 40)
    gt_total[i,"Term"] <- allRes[i, "Term"]
    gt_total[i,"NGenes"] <- dim(gt)[1]
    gt_genes[[allRes[i, "GO.ID"]]]<-gt
  }
  tab_gt <- gt_total[,2]
  names(tab_gt)<-gt_total[,1]
  
  main <- paste("frequency of ",GO," Terms",sep=" ")
  ylab <- paste(GO,"Terms",sep=" ")
  p <- barchart(tab_gt,main=main, ylab=ylab,  
                  xlab="Frequency",
                   panel=function(...){ 
                   args <- list(...)
                   panel.text(args$x, args$y, args$x, pos=4, offset=1)
                   panel.barchart(...)
                 }
)
  print(p)
  
  return(gt_genes)
}
