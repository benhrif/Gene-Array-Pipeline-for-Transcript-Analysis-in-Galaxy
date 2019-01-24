annoter_variants <- function(listeVariants){
  library("biomaRt")
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
  mart<-useDataset("hsapiens_gene_ensembl",mart)
  annot <- getBM(attributes=c("ensembl_transcript_id","ensembl_gene_id","hgnc_symbol","refseq_mrna","refseq_mrna_predicted"),
                 filters="ensembl_transcript_id",values=rownames(listeVariants), mart=mart)
  
  all <- merge(listeVariants,annot,by.x=0, by.y=1)
  return(all)
}
