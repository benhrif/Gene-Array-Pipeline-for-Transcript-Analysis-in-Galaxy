package.installed <- installed.packages()
suppressMessages(suppressWarnings(source("http://bioconductor.org/biocLite.R")))
if(!"VennDiagram" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("VennDiagram")))
}
if(!"methods" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("methods")))
}
if(!"S4Vectors" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("S4Vectors")))
}
if(!"siggenes" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("siggenes")))
}
if(!"RankProd" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("RankProd")))
}
if(!"RSQLite" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("RSQLite")))
}
if(!"DBI" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("DBI")))
}
if(!"proxy" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(install.packages("proxy")))
}
if(!"oligo" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("oligo")))
}
if(!"FactoMineR" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("FactoMineR")))
}
if(!"gplots" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("gplots")))
}
if(!"qvalue" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("qvalue")))
}
if(!"maanova" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("maanova")))
}
if(!"topGO" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("topGO")))
}
if(!"Rgraphviz" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("Rgraphviz")))
}
if(!"lattice" %in% rownames(package.installed)){   
  suppressMessages(suppressWarnings(biocLite("lattice")))
}
if(!"pvclust" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("pvclust")))
}
if(!"Biobase" %in% rownames(package.installed)){
  suppressMessages(suppressWarnings(biocLite("Biobase")))
}
