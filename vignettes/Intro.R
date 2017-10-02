## ------------------------------------------------------------------------
library(hiAnnotator)

## ------------------------------------------------------------------------
data(sites)
## sites object doesn't have a start & stop column to denote genomic range, hence soloStart parameter must be TRUE or a nasty error will be thrown!
alldata.rd <- makeGRanges(sites, soloStart = TRUE) 

data(genes)
## adding freeze populates SeqInfo slot of GRanges object.
genes.rd <- makeGRanges(genes, freeze = "hg18") 

## ---- eval = FALSE-------------------------------------------------------
#  refflat <- getUCSCtable("refFlat", "RefSeq Genes")
#  genes <- makeGRanges(refflat)

## ------------------------------------------------------------------------
nearestGenes <- getNearestFeature(alldata.rd, genes.rd, "NearestGene")
head(nearestGenes)

# nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene", parallel=TRUE)

## get nearest 5' genes
nearestGenes <- getNearestFeature(alldata.rd, genes.rd, "NearestGene", side = "5p") 
head(nearestGenes)

## get nearest 3' genes
nearestGenes <- getNearestFeature(alldata.rd, genes.rd, "NearestGene", side = "3p")
head(nearestGenes)

## get midpoint of genes
nearestGenes <- getNearestFeature(alldata.rd, genes.rd, "NearestGene", side = "midpoint")
head(nearestGenes)

### get two nearest upstream and downstream genes relative the query
nearestTwoGenes <- get2NearestFeature(alldata.rd, genes.rd, "NearestGene")
head(nearestTwoGenes)

## ------------------------------------------------------------------------
geneCounts <- getFeatureCounts(alldata.rd, genes.rd, "NumOfGene")
head(geneCounts)
# geneCounts <- getFeatureCounts(alldata.rd, genes.rd, "NumOfGene", parallel=TRUE)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  geneCounts <- getFeatureCounts(alldata.rd, genes.rd, "NumOfGene",
#                                 doInChunks = TRUE, chunkSize = 100)
#  head(geneCounts)
#  
#  geneCounts <- getFeatureCountsBig(alldata.rd, genes.rd, "NumOfGene")
#  head(geneCounts)

## ------------------------------------------------------------------------
## Shows which feature(s) a position was found in.
InGenes <- getSitesInFeature(alldata.rd, genes.rd, "InGene")
head(InGenes)

## Simply shows TRUE/FALSE 
InGenes <- getSitesInFeature(alldata.rd, genes.rd, "InGene", asBool = TRUE)
head(InGenes)

# InGenes <- getSitesInFeature(alldata.rd, genes.rd, "InGene", asBool=TRUE, parallel=TRUE)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  doAnnotation(annotType = "within", alldata.rd, genes.rd, "InGene")
#  doAnnotation(annotType = "counts", alldata.rd, genes.rd, "NumOfGene")
#  doAnnotation(annotType = "countsBig", alldata.rd, genes.rd, "ChipSeqCounts")
#  doAnnotation(annotType = "nearest", alldata.rd, genes.rd, "NearestGene")
#  doAnnotation(annotType = "twoNearest", alldata.rd, genes.rd, "TwoNearestGenes")
#  geneCheck <- function(x, wanted) { x$isWantedGene <- x$InGene %in% wanted;
#                                     return(x) }
#  doAnnotation(annotType = "within", alldata.rd, genes.rd, "InGene",
#               postProcessFun = geneCheck,
#               postProcessFunArgs = list("wanted" = c("FOXJ3", "SEPT9", "RPTOR")) )

## ---- eval=TRUE, echo=TRUE-----------------------------------------------
res <- doAnnotation(annotType = "within", alldata.rd, genes.rd, "InGene", asBool = TRUE)
plotdisFeature(res, "virus", "InGene")

res <- doAnnotation(annotType = "nearest", alldata.rd, genes.rd, "NearestGene", side = '5p')
plotdisFeature(res, "virus", "X5pNearestGeneDist")

data(sites.ctrl)
sites$type <- "expr"
sites <- rbind(sites,sites.ctrl)
alldata.rd <- makeGRanges(sites, soloStart = TRUE)
res <- doAnnotation(annotType = "within", alldata.rd, genes.rd, "InGene", asBool = TRUE)
plotdisFeature(res, "virus", "InGene")
plotdisFeature(res, "virus", "InGene", typeRatio = TRUE)

## ----par_examples, eval=FALSE, echo=TRUE---------------------------------
#  ## Example 1: library(doSMP)
#  w <- startWorkers(2)
#  registerDoSMP(w)
#  getNearestFeature(..., parallel = TRUE)
#  
#  ## Example 2: library(doMC)
#  registerDoMC(2)
#  getNearestFeature(..., parallel = TRUE)
#  
#  ## Example 3: library(doSNOW)
#  cl <- makeCluster(2, type = "SOCK")
#  registerDoSNOW(cl)
#  getNearestFeature(..., parallel = TRUE)
#  
#  ## Example 4: library(doParallel)
#  cl <- makeCluster(2)
#  registerDoParallel(cl)
#  getNearestFeature(..., parallel = TRUE)

