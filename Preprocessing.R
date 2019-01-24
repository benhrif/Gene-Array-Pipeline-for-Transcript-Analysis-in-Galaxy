
## load package
package.installed <- installed.packages()

suppressMessages(suppressWarnings(library("RSQLite")))
suppressMessages(suppressWarnings(library("DBI")))
suppressMessages(suppressWarnings(library("oligo")))

args <- suppressMessages(suppressWarnings(commandArgs(TRUE)));

suppressWarnings(suppressMessages(load(args[1])));
norm_method <- args[2]



####preprocessing

bg <- suppressMessages(suppressWarnings(backgroundCorrect(rawData)))
nm <- suppressMessages(suppressWarnings(normalize(bg, method = norm_method)))
ppData<-suppressMessages(suppressWarnings(rma(nm,background = FALSE, normalize = FALSE, target = "core")))

##output
suppressMessages(suppressWarnings(save(ppData,file="ppData.RData")))
