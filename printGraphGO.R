printGraphGO <- function(geneList,GO,Database,prefix){

affyLib<-Database

topDiffGenes <- function(geneList) {
  return(geneList)
}

sampleGOdata <- new("topGOdata",
                      description = "Simple session", ontology = GO,
                      allGenes = geneList,
                      geneSel = topDiffGenes,
                      annot = annFUN.db, affyLib = affyLib)

weight01.ks <- runTest(sampleGOdata, algorithm = "weight01", statistic = "ks")

test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test")

resultElim <- getSigGroups(sampleGOdata, test.stat)

allRes <- GenTable(sampleGOdata,  elim = resultElim,
                     weight01=weight01.ks,
                     orderBy = "elim",topNodes = 40)

printGraph(sampleGOdata, resultElim, firstSigNodes = 5,weight01.ks,
           fn.prefix = prefix, useInfo = "def", pdfSW = TRUE)

return(list(allRes,sampleGOdata))

}
