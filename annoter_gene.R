annoter_gene<-function(liste_gene,Database){
  chip <- substr(Database,1,nchar(Database)-3)
  SYMBOL <- get(paste(chip,"SYMBOL",sep=""))
  ACCNUM <- get(paste(chip,"ACCNUM",sep=""))
  GENENAME <- get(paste(chip,"GENENAME",sep=""))
  ENSEMBL <- get(paste(chip,"ENSEMBL",sep=""))
  suppressMessages(library(Database,character.only = TRUE))
  
  my_frame<-data.frame(liste_gene)
  Annot <- data.frame(SYMBOL=sapply(contents(SYMBOL),
                            paste, collapse=", "),
                      ACCNUM=sapply(contents(ACCNUM),
                                     paste, collapse=", "), 
                       DESC=sapply(contents(GENENAME), 
                                   paste, collapse=", "),
                       ENSEMBL=sapply(contents(ENSEMBL), 
                                     paste, collapse=", "))
  all <- merge(my_frame,Annot,by.x=0, by.y=0)
  all <- all[which(all[,"SYMBOL"]!="NA"),]
  return(all)
}
