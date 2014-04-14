#' Annotating RangedData or GRanges objects with hiAnnotator.
#'
#' hiAnnotator contains set of functions which allow users to annotate a RangedData or GRanges object with custom set of annotations. The basic philosophy of this package is to take two RangedData or GRanges objects (query & subject) with common set of space/seqnames (i.e. chromosomes) and return associated annotation per space/seqnames and rows from the query matching space/seqnames and rows from the subject (i.e. genes or cpg islands).
# The package comes with three types of annotation functions which calculates if a position from query is: within a feature, near a feature, or count features in defined window sizes. Moreover, one can utilize parallel backend for each annotation function to utilize the foreach package. In addition, the package is equipped with wrapper functions, which finds appropriate columns needed to make a RangedData or GRanges object from a common dataframe..
#'
#' @import GenomicRanges foreach iterators rtracklayer plyr BSgenome
#' @docType package
#' @name hiAnnotator
NULL

#' Sample Retrovirus Integration Sites data
#' 
#' A sample dataset containing collection of unique HIV & MLV integration sites in the human genome mapped to UCSC freeze hg18 from PMID: 12805549.
#' 
#' \itemize{
#'   \item Sequence. Name of the DNA sequence which was aligned to the host genome. This is also a unique ID.
#'   \item Position. The genomic coordinate of the integration site.
#'   \item Chr. The chromosome of the integration site. 
#'   \item Ort. The orientation or strand of the integration site. 
#'   \item virus. Name of the virus used for the experiment and a given sequencing clone.
#' }
#' 
#' @source \url{http://www.ncbi.nlm.nih.gov/pubmed/?term=12805549}
#' @docType data
#' @keywords datasets
#' @format A data frame with 1303 rows and 5 variables
#' @name sites
NULL

#' Sample RefSeq genes annotation
#' 
#' A sample annotation containing collection of genes from RefSeq database in the human genome mapped to UCSC freeze hg18. See UCSC table description page for the details regarding the column headings.
#' 
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgTables?db=hg18&hgta_table=refGene&hgta_doSchema=describe+table+schema}
#' @docType data
#' @keywords datasets
#' @format A data frame with 33965 rows and 9 variables
#' @name genes
NULL

#' Initiate UCSC genome browser session given the freeze argument.
#'
#' @param freeze one of following: hg18, mm8, rheM, etc. Default is hg18.
#'
#' @return browser session object compatible with rtracklayer functions.
#'
#' @seealso \code{\link{getUCSCtable}}, \code{\link{makeRangedData}}, 
#' \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, 
#' \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' #session <- makeUCSCsession()
#' #genome(session)
#' #session <- makeUCSCsession("mm8")
#' #genome(session)
makeUCSCsession <- function(freeze="hg18") {
  bsession <- browserSession()
  genome(bsession) <- freeze
  bsession
}

#' Obtain a UCSC annotation table given the table & track name.
#'
#' @param tableName Name of the annotation table as it appears on UCSC browser.
#' @param trackName Name of the track annotation table appears in on UCSC browser.
#' @param bsession UCSC session object returned by \code{\link{makeUCSCsession}} or \code{\link{browserSession}}. If left NULL the function will call \code{\link{makeUCSCsession}} with the provided freeze to initiate a session.
#' @param freeze one of following: hg18, mm8, rheM, etc. Default is hg18.
#' @param ... Arguments to be passed to \code{\link{ucscTableQuery}}.
#'
#' @return a dataframe containing the annotation data.
#'
#' @seealso \code{\link{makeUCSCsession}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' #refflat <- getUCSCtable("refFlat","RefSeq Genes") ## same as session <- makeUCSCsession(); 
#' #refflat <- getUCSCtable("refFlat", "RefSeq Genes", bsession=session, freeze="hg18")
#'
getUCSCtable <- function(tableName, trackName, 
                         bsession=NULL, freeze="hg18", ...) {
  if(is.null(bsession)) { 
    bsession <- makeUCSCsession(freeze)
  }
  if(!tableName %in% tableNames(ucscTableQuery(bsession,track=trackName))) {
    stop(paste("The provided table name:",tableName,"doesn't exists in track",
               trackName,"on UCSC for",freeze,"genome"))
  }
  
  ## using getTable() instead of track() due to "No supported output types" error 
  ## for certain annotation types.    
  getTable(ucscTableQuery(bsession,track=trackName,table=tableName,...))
}

#' Find the column index of interest given the potential choices.
#'
#' The function finds relevant column(s) of interest from a vector of column names derived from a dataframe or a RangedData object. If no usable column is found, the function spits out a relevant error or returns the index of the usable column(s). This is an assistant function called by functions listed in the see also section.
#'
#' @param col.names column names from a dataframe or a RangedData object
#' @param col.options potential column names or partial names that may exist in col.names
#' @param col.type type of column information the function is searching for, used in construction of error messages. Default is NULL.
#' @param multiple.ok if multiple matches are found then return indices, else spit an error out. Default is TRUE.
#'
#' @return the index of usable column(s) or an error if no applicable column is found.
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' data(sites)
#' names(sites)
#' getRelevantCol(names(sites),c("chr","chromosome","tname","space","chrom","contig"),"space")
#' getRelevantCol(names(sites),c("ort","orientation","strand"),"strand")
getRelevantCol <- function(col.names, col.options, 
                           col.type=NULL, multiple.ok=FALSE) {
  answer <- unique(as.numeric(unlist(sapply(col.options,
                                            function(x) grep(x,col.names,
                                                             ignore.case=TRUE))
                                     )))
  if(length(answer)>1) {
    if(!multiple.ok) {
      stop(paste("More than one",col.type,"based column found:",
                 paste(col.names[answer],sep="",collapse=", ")))
    } else {
      answer
    }
  } else if (length(answer)==0) {
    stop(paste("No",col.type,"based column found."))
  } else {
    answer
  }
}

#' Breaks two GRanges/RangedData objects into chunks of N size.
#'
#' Given a query and subject GRanges/RangedData objects, the function breaks query into chunks of N size where each chunk has a respective subject object filtered by seqnames/space present in the query chunk. This is a helper function used by one of the annotation function in 'See Also' section where each chunk is sent to a parallel node for processing.
#'
#' @param sites.rd a RangedData/GRanges object.
#' @param features.rd a RangedData/GRanges object.
#' @param chunkSize number of rows to use per chunk of query. Default to length(sites.rd)/detectCores() or length(query)/getDoParWorkers() depending on parallel backend registered.
#'
#' @return a list of RangedData/GRanges objects where each element is of length 2 representing query & subject chunks. 
#'
#' @seealso \code{\link{makeGRanges}}, \code{\link{doAnnotation}},
#' \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}, 
#' \code{\link{getFeatureCounts}}.
#'
#' @export
#'
#' @examples
#' data(sites)
#' data(genes)
#' sites <- makeGRanges(sites,soloStart=TRUE)
#' genes <- makeGRanges(genes)
#' makeChunks(sites, genes)
makeChunks <- function(sites.rd, features.rd, chunkSize=NULL) {
  # do a quick check of things
  .checkArgsSetDefaults()
  rm("query","subject")
  
  # make chunks
  chunks <- breakInChunks(length(sites.rd), 
                          ifelse(!is.null(chunkSize), 
                                 length(sites.rd)/chunkSize,
                                 ifelse(!is.null(is.null(getDoParWorkers())),
                                        length(sites.rd)/getDoParWorkers(),
                                        length(sites.rd)/detectCores())))
  
  mapply(function(x,y) {
    new.query <- sites.rd[x:y,]
    new.query <- keepSeqlevels(new.query, 
                               value=unique(as.character(seqnames(new.query))))
    new.subject <- suppressWarnings(keepSeqlevels(features.rd, 
                                                  value=seqlevels(new.query)))
    if(RangedDataFlag) {
      new.query <- as(new.query,"RangedData")
      new.subject <- as(new.subject,"RangedData")
    }
    list("query"=new.query, "subject"=new.subject)
  }, start(chunks), end(chunks), SIMPLIFY=FALSE, USE.NAMES=FALSE)
}

#' Clean the supplied string from punctuations and spaces.
#'
#' Function to clean the supplied string from punctuations and spaces so it can be used as column headings.
#'
#' @param x string or a vector to be cleaned.
#' @param description OPTIONAL string identifying the purpose of the supplied string in x to be displayed in the cleaning message. This triggers a message.
#'
#' @return cleaned string or a vector.
#'
#' @seealso \code{\link{getFeatureCounts}}, \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' cleanColname("HIV-test")
#' cleanColname("HIV*test")
#' cleanColname("HIV-test","myAlias")
cleanColname <- function(x, description=NULL) {
  newname <- gsub("[._]+","_",make.names(x,unique=TRUE))
  if(any(newname != x ))
    if(!is.null(description))
      message("Cleaning the supplied '",description,"'")
  newname
}

#' Resize a RangedData object.
#'
#' Function to resize a RangedData object by the given width, max space/chromosome size, and the boundary. This is one of the helper function used in \code{\link{getFeatureCountsBig}}. 
#'
#' @param rd a RangedData object
#' @param width the width of the resized ranges.
#' @param boundary same as fix parameter in \code{\link[IRanges]{resize}}. One of "start", "end", and "center". Default is "center".
#' @param spaceSizes named vector of chromosome/space sizes to be used for fixing ranges off the limits.
#' @param spaceMin lowest value to be allowed as the start position. Default is 1.
#' @param limitLess whether to limit the resized ranges by spaceSizes and spaceMin. Default is FALSE. In genomic context, this means starts can be <=0 and ends can be > chromosome size depending on the width. 
#'
#' @return a RangedData object with ranges sized to width and/or truncated to spaceSizes and spaceMin unless limitLess=TRUE.
#'
#' @seealso \code{\link{getFeatureCountsBig}}, \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#' seqLengths <- structure(c(247249719L, 242951149L, 199501827L, 191273063L, 180857866L,170899992L, 158821424L, 146274826L, 140273252L, 135374737L, 134452384L, 132349534L, 114142980L, 106368585L, 100338915L, 88827254L, 78774742L, 76117153L, 63811651L, 62435964L, 46944323L, 49691432L, 154913754L, 57772954L, 16571L), .Names = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"))
#' resizeRangedData(alldata.rd,width=10000,spaceSizes=seqLengths)
#' resizeRangedData(alldata.rd,width=10000,limitLess=TRUE)
resizeRangedData <- function(rd, width=NULL, 
                             boundary="center", spaceSizes=NULL, 
                             spaceMin=1, limitLess=FALSE) {
  stopifnot(!is.null(width))
  
  RangedData <- FALSE
  if(is(rd,"RangedData")) {
    RangedData <- TRUE
  }
  
  new.rd <- resize(ranges(rd),width,fix=boundary)
  new.rd <- restrict(new.rd, start=spaceMin)
  if(!as.logical(limitLess)) {
    stopifnot(!is.null(spaceSizes))    
    overEdge <- sapply(names(new.rd),
                       function(x) which(end(new.rd[[x]])>spaceSizes[x]),
                       simplify=FALSE)
    if (any(sapply(overEdge,length)>0)) {
      overEdge.chr <- names(which(sapply(overEdge,length)>0))
      for (f in overEdge.chr) {
        end(new.rd[[f]])[ as.numeric(unlist(overEdge[f]))] <- spaceSizes[f]
      } 
    }  	
  }
  
  if(RangedData) {
    new.rd <- as(new.rd, "RangedData")
    values(new.rd) <- values(rd)
  } else {
    ranges(rd) <- new.rd
    new.rd <- rd
  }
  
  return(new.rd)
}

#' Make a sorted RangedData object from a dataframe. 
#'
#' The function converts a dataframe into a RangedData object without too much hassle of renaming column names. The function finds column names that sound like space, chromosome, start, stop, position, etc and puts them in respective slots to facilitate the conversion of a dataframe to a RangedData object. If more than one column that sounds like start, stop, or position is present, the function will use the first match as the representative. It is recommended to run this function before utilizing any other annotation functions since it will sort the object by chromosome and position for copying annotations back to their respective rows confidently. It is recommended to use a \code{\link[GenomicRanges]{GRanges}} object instead of a \code{\link[IRanges]{RangedData}} object if number of distinct chromosomes or targets is greater than 50. 
#'
#' @param x dataframe to be converted into a RangedData object
#' @param positionsOnly boolean flag indicating to return only position based data or everything from the dataframe. Defaults to FALSE.
#' @param soloStart flag denoting whether only one position based column is available. In other words, only starts are present and no stops. Default=FALSE.
#' @param chromCol use the defined column name for space/chromosome based data from the dataframe. Defaults to NULL.
#' @param strandCol use the defined column name for strand or orientation from the dataframe. Defaults to NULL.
#' @param startCol use the defined column name for start coordinate from the dataframe. Defaults to NULL.
#' @param stopCol use the defined column name for stop coordinate from the dataframe. Defaults to NULL and not required if soloStart=TRUE.
#' @param asGRanges make a GRanges object instead. This is useful if number of target or chromosomes is greater than 50. Default is FALSE.
#'
#' @return a RangedData object converted from x.
#'
#' @seealso \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData object
#' data(sites)
#' head(sites)
#'
#' makeGRanges(sites,soloStart=TRUE)
#' makeRangedData(sites,soloStart=TRUE,positionsOnly=TRUE)
#' # makeRangedData(sites) # should yield an error
#'
#' #data(genes)
#' #head(genes)
#'
#' #makeRangedData(genes,soloStart=TRUE)
#' #makeGRanges(genes)
makeRangedData <- function(x, positionsOnly=FALSE, soloStart=FALSE, chromCol=NULL, 
                           strandCol=NULL, startCol=NULL, stopCol=NULL, 
                           asGRanges=FALSE) {
  ## set column names for space and strand if not provided ##
  if(is.null(chromCol)) {
    colIndex <- getRelevantCol(names(x),
                               c("chr","chromosome","tname","space","chrom",
                                 "contig","seqnames"),
                               "space")
    chromCol <- names(x)[colIndex]
  }
  x$space <- x[,chromCol]
  
  if(is.null(strandCol)) {
    colIndex <- getRelevantCol(names(x),
                               c("ort","orientation","strand"),
                               "strand")
    strandCol <- names(x)[colIndex]
  }    
  x$strand <- x[,strandCol]
  
  if(is.null(startCol)) {
    startCol <- getRelevantCol(names(x),
                               c("position", "intsite", "txstart",
                                 "start", "chromstart"),
                               "start", multiple.ok=TRUE)
    startCol <- names(x)[startCol[1]]
  }
  x$start <- x[,startCol]
  
  ## only do stop if soloStart=F ##
  if(!as.logical(soloStart) & is.null(stopCol)) {
    stopCol <- getRelevantCol(names(x),
                              c("txend", "end", "stop", "chromend"),
                              "end", multiple.ok=TRUE)
    stopCol <- names(x)[stopCol[1]]
  }
  x$end <- x[,stopCol]
  
  ## do some testing for NAs in space, start or end ##
  if(any(is.na(x$start))) {
    stop("NAs found in column containing start positions")
  }
  
  if(any(is.na(x$space))) {
    stop("NAs found in column containing chromosome or space information")
  }
  
  ## convert any factor columns to character to avoid downstream issues ##
  factorCols <- sapply(x,is.factor)
  if(any(factorCols)) {
    for (y in names(which(factorCols))) { 
      x[,y] <- as.character(x[,y])
      if(!any(is.na(suppressWarnings(as.numeric(x[,y]))))) { 
        x[,y] <- as.numeric(x[,y]) 
      }
    }        
  }
  
  ## if start and end coordinates are present, sort by midpoint ##
  ## else if only single coordinate is present, then add the end column and sort ##
  if(length(startCol)>0 & length(stopCol)>0) {       
    if(any(is.na(x$end))) {
      stop("NAs found in column containing end positions")
    }
    x$mid <- with(x,(start+end)/2)
    x <- arrange(x,space,mid)
    x$mid <- NULL
  } else {  
    x <- arrange(x,space,start)
    x$end <- x$start
  }
  
  if(as.logical(positionsOnly)) {
    x <- x[,c("space","start","end","strand")]    
  } 
  
  metadataCols <- grep("space|start|end|strand", names(x),
                       invert=TRUE, value=TRUE, fixed=FALSE)
  metadataCols <- metadataCols[!is.na(metadataCols)]
  if(asGRanges) {
    sites.rd <- GRanges(seqnames=x$space, IRanges(start=x$start, end=x$end),
                        strand=x$strand)
    ## Loop through incase only one metadataCol is present which returns a vector 
    ## instead of a dataframe...DataFrame(x[,metadataCols]) may not work! 
    for(f in metadataCols) {       
      mcols(sites.rd)[[f]] <- x[,f]
    }
  } else {
    sites.rd <- RangedData(space=x$space, IRanges(start=x$start, end=x$end),
                        strand=x$strand)
    for(f in metadataCols) {
      sites.rd[[f]] <- x[,f]
    }
  }
  
  sites.rd
}

#' Make a sorted GRanges object from a dataframe. 
#'
#' The function converts a dataframe into a GRanges object without too much hassle of renaming column names. The function finds column names that sound like seqname, chromosome, start, stop, position, etc and puts them in respective slots to facilitate the conversion of a dataframe to a GRanges object. If more than one column that sounds like start, stop, or position is present, the function will use the first match as the representative. It is recommended to run this function before utilizing any other annotation functions since it will sort the object by chromosome and position for copying annotations back to their respective rows confidently. This function wraps around \code{\link{makeRangedData}} with asGRanges parameter enabled. It is recommended to use a \code{\link[GenomicRanges]{GRanges}} object instead of a \code{\link[IRanges]{RangedData}} object if number of distinct chromosomes or targets is greater than 50. 
#'
#' @param x dataframe to be converted into a GRanges object
#' @param freeze UCSC genome version of the data in x. Default is NULL. This parameter is generally used to populate \code{\link{seqinfo}} slot of GRanges objects.
#' @param ... parameters for \code{\link{makeRangedData}} except for 'asGRanges'.
#'
#' @return a GRanges object converted from x.
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{getNearestFeature}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(genes)
#'
#' makeGRanges(genes, soloStart=TRUE)
#' makeGRanges(genes)
#' #makeGRanges(genes, freeze="hg18", soloStart=TRUE)
#' #makeGRanges(genes, freeze="hg18")
makeGRanges <- function(x, freeze=NULL, ...) {    
  sites.gr <- makeRangedData(x, asGRanges=TRUE, ...)
  
  if(!is.null(freeze)) {
    genomeLib <- grep(freeze,installed.genomes(),value=TRUE)
    if(length(genomeLib)!=0) {		
      bsGenomeObject <- strsplit(genomeLib,"\\.")[[1]][2]
      chrom.info <- seqlengths(do.call(`:::`,list(genomeLib,bsGenomeObject)))
    } else {
      ## get the chromInfo file from UCSC
      z <- gzcon(url(paste0("http://hgdownload.cse.ucsc.edu/goldenPath/", 
                            freeze, "/database/chromInfo.txt.gz")))
      zlines <- try(readLines(z))
      close(z)
      if (class(zlines)=="try-error") 
        stop("Could not get thru to UCSC server - try later or drop the freeze parameter!")
      raw.data <- textConnection(zlines)
      chrom.info <- read.delim(raw.data, header=FALSE,
                               stringsAsFactors=FALSE)[,1:2]
      chrom.info <- structure(chrom.info$V2, names=chrom.info$V1)
      close(raw.data)
    }
    
    # append to sites.gr as Seqinfo slot & reorder things while at it!
    chrom.info <- chrom.info[names(chrom.info) %in% levels(seqnames(sites.gr))]
    chrom.info <- chrom.info[order(suppressWarnings(
      as.numeric(gsub("chr","",names(chrom.info))))
    )]
    seqinfo(sites.gr,
            new2old=order(suppressWarnings(
              as.numeric(gsub("chr","",levels(seqnames(sites.gr))))
            ))) <- 
      Seqinfo(names(chrom.info), chrom.info, 
              rep(FALSE,length(chrom.info)), freeze)
  }
  
  sites.gr
}

#' Get nearest annotation boundary for a position range. 
#'
#' Given a query object, the function retrieves the nearest feature and its properties from a subject and then appends them as new columns within the query object. When used in genomic context, the function can be used to retrieve a nearest gene 5' or 3' end relative to genomic position of interest.
#'
#' @param sites.rd RangedData/GRanges object to be used as the query.
#' @param features.rd RangedData/GRanges object to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves a core!
#' @param side boundary of annotation to use to calculate the nearest distance. Options are '5p','3p', 'either'(default), or 'midpoint'.
#' @param feature.colnam column name from features.rd to be used for retrieving the nearest feature name. By default this is NULL assuming that features.rd has a column that includes the word 'name' somewhere in it.
#' @param dists.only flag to return distances only. If this is TRUE, then 'feature.colnam' is not required and only distance to the nearest feature will be returned. By default this is FALSE.
#' @param parallel use parallel backend to perform calculation with \code{\link[foreach]{foreach}}. Defaults to FALSE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link[foreach]{registerDoSEQ}}.
#' @param relativeTo calculate distance relative to query or subject. Default is 'subject'. This essentially means whether to use query or subject as the anchor point to get distance from!
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd.
#'
#' @note
#' \itemize{
#'   \item When side='midpoint', the distance to nearest feature is calculated by (start+stop)/2. 
#'   \item Try not to use this function for >50 spaces/seqnames/chromosomes unless you have tons fo memory. 
#'   \item If strand information doesn't exist, then everything is defaulted to '+' orientation (5' -> 3')
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered before running the function. One can use any of the following libraries compatible with \code{\link[foreach]{foreach}}: doMC, doSMP, doSNOW, doMPI, doParallel. For example: library(doMC); registerDoMC(2)
#'   \item When relativeTo="subject", the biological distance is relative to subject, meaning, the function reports the distance to query from subject (i.e. an integration site is upstream or downstream from a gene). When relativeTo="query", the distance is from the point of view of query or an integration site (i.e. gene is upstream or downstream from an integration site).
#' }
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}, \code{\link{get2NearestFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene")
#' nearestGenes
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",side="5p")
#' nearestGenes
#' \dontrun{
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",side="3p")
#' nearestGenes
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",side="midpoint")
#' ## Parallel version of getNearestFeature
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene", parallel=TRUE)
#' nearestGenes
#' }
getNearestFeature <- function(sites.rd, features.rd, 
                              colnam=NULL, side="either", feature.colnam=NULL, 
                              dists.only=FALSE, parallel=FALSE, 
                              relativeTo='subject') {
  
  .checkArgsSetDefaults()

  if(!dists.only) {    
    mcols(subject)$featureName <- mcols(features.rd)[,feature.colnam]
  }
  rm(features.rd)
  
  if(side %in% c('5p','3p','midpoint')) {
    ##get only 5 prime sides of features
    if (side=='5p')
      subject <- flank(subject, width=-1) 
    
    ##get only 3 prime sides of features
    if (side=='3p')
      subject <- flank(subject, width=-1, start=FALSE) 
    
    ##get (start+stop)/2 of features
    if (side=='midpoint')
      ranges(subject) <- IRanges(mid(ranges(subject)), width=1)
  }

  prefix <- ifelse(side=="either","",side)    
  colnam <- cleanColname(colnam)
  
  ## chunksize the objects for parallel processing ##
  chunks <- if(parallel) { 
    makeChunks(query, subject)
  } else {
    list(list("query"=query, "subject"=subject))
  }

  ## first get the nearest indices, respective tempyIDs, and distances ##  
  res <- foreach(x=iter(chunks), .inorder=FALSE, 
                 .export=c("side","relativeTo"),
                 .packages=c("GenomicRanges","plyr")) %dopar% { 
                   res.x <- as.data.frame(nearest(x$query, x$subject, 
                                                  select="all", 
                                                  ignore.strand=TRUE))
                   res.x$qID <- mcols(x$query)$tempyID[res.x$queryHits]
                   res.x$sID <- mcols(x$subject)$tempyID[res.x$subjectHits]
                   res.x <- getLowestDists(x$query, x$subject, res.x, 
                                           side, relativeTo)
                   counts <- count(res.x,"queryHits")
                   merge(res.x, counts)
                 }

  if(!dists.only) {    
    ## for the feature of shortest indices, get the names, and strand attributes    
    ## fix cases where >1 equally nearest features were returned by concatenating 
    ## feature names and strand while returning one distance per query
    
    res <- foreach(x=iter(res), y=iter(sapply(chunks,"[[","subject")), 
                   .inorder=FALSE, .combine=rbind) %dopar% {
                     # make sure x & y have the respective data chunks! #
                     stopifnot(all(x$sID %in% mcols(y)$tempyID))
                     
                     x$featureName <- mcols(y)[,"featureName"][x$subjectHits]
                     x$strand <- as.character(strand(y))[x$subjectHits]
                     
                     # isolate non-singletons to save time & memory! #
                     besties <- droplevels(subset(x,freq==1))
                     x <- droplevels(subset(x,freq>1))
                     
                     # use tapply instead of ddply() or by() because it's a lot faster on larger datasets #
                     bore <- with(x, 
                                  sapply(tapply(featureName,queryHits,unique), 
                                         paste, collapse=","))
                     x$featureName <- as.character(bore[as.character(x$queryHits)])
                     
                     bore <- with(x, sapply(tapply(strand,queryHits,unique), 
                                            paste, collapse=","))
                     x$strand <- as.character(bore[as.character(x$queryHits)])
                     
                     x <- unique(x[,c("queryHits", "qID", "dist",
                                      "featureName","strand")])
                     
                     # put singletons & curated non-singletons back together! #
                     besties <- rbind(besties[,names(x)], x)
                     besties <- arrange(besties,qID)                                          
                     
                     besties
                   }
    ## change column names for swift merging by .mergeAndReturn() ##
    names(res)[grepl("featureName",names(res))] <- paste0(prefix,colnam)
    names(res)[grepl("strand",names(res))] <- paste0(prefix,colnam,"Ort")

    } else {
    ## fix cases where >1 equally nearest features were returned by choosing 1 distance
    res <- foreach(x=iter(res), .inorder=FALSE, 
                   .combine=rbind ) %dopar% {
                     unique(x[,c("queryHits","qID","dist")])                     
                   }    
  }
  
  rm(chunks)

  ## change distance column name for .mergeAndReturn() ##
  names(res)[grepl("dist",names(res))] <- paste0(prefix,colnam,"Dist")
  
  # Do a last check to make sure there is only 1 hit per qID #
  # This is useful in cases where two equally nearest distances 
  # but in opposite directions are returned #                     
  test <- duplicated(res$qID)
  if(any(test)) {                       
    res <- res[!test,]
  }
  
  ## merge results to the query object and return it ##
  .mergeAndReturn()  
  
  sites.rd
}

#' Get two nearest upstream and downstream annotation boundary for a position range. 
#'
#' Given a query object, the function retrieves the two nearest feature upstream and downstream along with their properties from a subject and then appends them as new columns within the query object. When used in genomic context, the function can be used to retrieve two nearest gene upstream and downstream of the genomic position of interest.
#'
#' @param sites.rd RangedData/GRanges object to be used as the query.
#' @param features.rd RangedData/GRanges object to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves a core!
#' @param side boundary of annotation to use to calculate the nearest distance. Options are '5p','3p', 'either'(default), or 'midpoint'.
#' @param feature.colnam column name from features.rd to be used for retrieving the nearest feature name. By default this is NULL assuming that features.rd has a column that includes the word 'name' somewhere in it.
#' @param relativeTo calculate distance relative to query or subject. Default is 'subject'. See documentation of  \code{\link{getNearestFeature}} for more information.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd.
#'
#' @note
#' \itemize{
#'   \item When side='midpoint', the distance to nearest feature is calculated by (start+stop)/2. 
#'   \item For cases where a position is at the edge and there are no feature up/down stream since it would fall off the chromosome, the function simply returns 'NA'. 
#'   \item If there are multiple locations where a query falls into, the function arbitrarily chooses one to serve as the nearest feature, then reports 2 upstream & downstream feature. That may occasionally yield features which are the same upstream and downstream, which is commonly encountered when studying spliced genes or phenomena related to it. 
#'   \item Try not to use this function for >50 spaces/seqnames/chromosomes unless you have tons fo memory. 
#'   \item If strand information doesn't exist, then everything is defaults to '+' orientation (5' -> 3')
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered before running the function. One can use any of the following libraries compatible with \code{\link[foreach]{foreach}}: doMC, doSMP, doSNOW, doMPI, doParallel. For example: library(doMC); registerDoMC(2)
#' }
#'
#' @seealso \code{\link{getNearestFeature}}, \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene")
#' nearestGenes
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene",side="5p")
#' nearestGenes
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene",side="3p")
#' nearestGenes
#' }
get2NearestFeature <- function(sites.rd, features.rd, 
                               colnam=NULL, side="either", 
                               feature.colnam=NULL, relativeTo="subject") {
  
  .checkArgsSetDefaults()
   
  ## make sure features.rd/subject is sorted ##
  mcols(subject)$featureName <- mcols(features.rd)[,feature.colnam]
  subject <- sort(subject)
  rm(features.rd)
  
  if(side %in% c('5p','3p','midpoint')) {
    ##get only 5 prime sides of features
    if (side=='5p')
      subject <- flank(subject, width=-1) 
    
    ##get only 3 prime sides of features
    if (side=='3p')
      subject <- flank(subject, width=-1, start=FALSE) 
    
    ##get (start+stop)/2 of features
    if (side=='midpoint')
      ranges(subject) <- IRanges(mid(ranges(subject)), width=1)
  }
  
  ## u = upstream, d = downstream
  ## thinking concept: u2.....u1.....intSite(+).....d1.....d2
  ## thinking concept: d2.....d1.....intSite(-).....u1.....u2
  ## searching concept: res.left2.....res.left1.....res....intSite....res.....res.right1.....res.right2
  
  ## first get the nearest indices, respective tempyIDs ##  
  res <- as.data.frame(nearest(query, subject, select="all", 
                               ignore.strand=TRUE))
  res$qID <- mcols(query)$tempyID[res$queryHits]
  res$qStrand <- as.character(strand(query))[res$queryHits]
  res <- getLowestDists(query, subject, res, side, relativeTo)
  
  ## perform upstream-downstream checks by testing distances  
  res$u2 <- with(res, ifelse(dist<0, 
                             ifelse(qStrand=="+", subjectHits - 1, 
                                    subjectHits + 1), 
                             ifelse(qStrand=="+", subjectHits - 2, 
                                    subjectHits + 2)))

  res$u1 <- with(res, ifelse(dist<0, subjectHits, 
                             ifelse(qStrand=="+", subjectHits - 1, 
                                    subjectHits + 1)))
  
  res$d1 <- with(res, ifelse(dist<0, 
                             ifelse(qStrand=="+", subjectHits + 1, 
                                    subjectHits - 1), 
                             subjectHits))
  
  res$d2 <- with(res, ifelse(dist<0, 
                             ifelse(qStrand=="+", subjectHits + 2, 
                                    subjectHits - 2), 
                             ifelse(qStrand=="+", subjectHits + 1, 
                                    subjectHits - 1)))
  
  prefix <- ifelse(side=="either","Either",side)
  
  message("u = upstream, d = downstream")
  message("thinking concept: u2.....u1.....intSite(+).....d1.....d2")
  message("thinking concept: d2.....d1.....intSite(-).....u1.....u2")
  
  colnam <- cleanColname(colnam)
  
  ## add columns back to query object
  all.res <- lapply(c("u1","u2","d1","d2"), function(f) {
    message(f)
    res.nrst <- res[,c("queryHits","subjectHits","qID",f)]
    
    # make sure we haven't jumped a chromosome by shifting nearest indices above #
    # if we did, then record it as NA at later a stage #    
    
    # fix cases where chosen subjectHits are off the length of subject object #
    fixed <- with(res.nrst, ifelse(get(f) < 1 | get(f) > length(subject), 
                                   subjectHits, get(f)))
    
    # do the chromosome test & tag rows which were off the subject length #
    res.nrst$qChr <- as.character(seqnames(query))[res.nrst$queryHits]
    res.nrst$sChr <- as.character(seqnames(subject))[fixed]
    rows <- res.nrst$qChr!=res.nrst$sChr | res.nrst[,f] < 1 | 
      res.nrst[,f] > length(subject)
    
    # overwrite subjectHits indices with that of interested motif for later steps #
    res.nrst$subjectHits <- res.nrst[,f]    
    
    # remove unnecessary columns #
    res.nrst[,f] <- NULL
    res.nrst$qChr <- NULL
    res.nrst$sChr <- NULL
    
    # extract cases which fell off the chromosome but drop any queryHits which found
    # multiple nearest hits and only one of them happened to be off the chromosome! 
    res.nrst.bad <- droplevels(res.nrst[rows & !res.nrst$queryHits %in% 
                                          res.nrst$queryHits[!rows],])
    res.nrst <- droplevels(res.nrst[!rows,])
    
    res.nrst <- getLowestDists(query, subject, res.nrst, side, relativeTo)
    res.nrst$featureName <- mcols(subject)[,"featureName"][res.nrst$subjectHits]
    res.nrst$strand <- as.character(strand(subject))[res.nrst$subjectHits]    
      
    res.nrst <- with(res.nrst,
                  by(res.nrst, queryHits, function(n) 
                    with(n, 
                         data.frame(queryHits=unique(queryHits), 
                                    qID=unique(qID), 
                                    dist=unique(dist), 
                                    featureName=paste(unique(featureName),
                                                      collapse=","),
                                    strand=paste(unique(strand),collapse=","),
                                    stringsAsFactors=FALSE)))
    )                     
    res.nrst <- do.call(rbind, res.nrst)
    
    # add back rows which fell off the edge of chromosome #
    if(any(rows)) {
      res.nrst.bad[,c("dist", "featureName", "strand")] <- NA
      res.nrst <- rbind(res.nrst, unique(res.nrst.bad[,names(res.nrst)]))
      res.nrst <- arrange(res.nrst, qID)
    }
    
    if (f == "u1") { coldef <- paste(prefix,colnam,"upStream1",sep=".") }
    if (f == "u2") { coldef <- paste(prefix,colnam,"upStream2",sep=".") }
    if (f == "d1") { coldef <- paste(prefix,colnam,"downStream1",sep=".") }
    if (f == "d2") { coldef <- paste(prefix,colnam,"downStream2",sep=".") }
    
    ## add meta columns to the result ##
    names(res.nrst)[grepl("featureName",names(res.nrst))] <- coldef
    names(res.nrst)[grepl("strand",names(res.nrst))] <- paste(coldef,"Ort",
                                                              sep=".")
    names(res.nrst)[grepl("dist",names(res.nrst))] <- paste(coldef,"Dist",
                                                            sep=".")
    
    res.nrst
  })
  
  res <- do.call(cbind, all.res)

  ## merge results to the query object and return it ##
  .mergeAndReturn()  
  
  sites.rd
}

#' Get the lowest biological distance from the 5' or 3' boundaries of query and subject.
#'
#' Given a query and subject with indicies from \code{\link[IRanges]{nearest}}, calculate the shortest biological distance to either boundaries of the query and subject. This is a helper function utilized in \code{\link{getNearestFeature}}, \code{\link{get2NearestFeature}}
#'
#' @param query GRanges object to be used as the query which holds data for 'queryHits' attribute of res.nrst.
#' @param subject GRanges object to be used as the subject which holds data for 'subjectHits' attribute of res.nrst.
#' @param res.nrst a dataframe of nearest indices as returned by \code{\link[IRanges]{nearest}}.
#' @param side boundary of subject/annotation to use to calculate the nearest distance. Options are '5p','3p', or the default 'either'.
#' @param relativeTo calculate distance relative to query or subject. Default is 'subject'. See documentation of  \code{\link{getNearestFeature}} for more information.
#'
#' @return res.nrst with lowest distances appended at the end.
#'
#' @note for cases where a query has multiple nearest neighbors or overlaps with >1 subjects, the function will choose the subject with the lowest absolute distance.
#'
#' @seealso \code{\link{getNearestFeature}}, \code{\link{get2NearestFeature}}.
#'
#' @export
getLowestDists <- function(query=NULL, subject=NULL, res.nrst=NULL, 
                           side="either", relativeTo="subject") {
  if(is.null(query) | is.null(subject) | is.null(res.nrst)) {
    stop("One of following is null: query, subject, res.nrst")
  }    
    
  if (side=="either") {
    ## get the lowest dist to either annot boundary from 5p side of the query
    dist.s <- start(query)[res.nrst$queryHits] - 
      start(subject)[res.nrst$subjectHits]
    dist.e <- start(query)[res.nrst$queryHits] - 
      end(subject)[res.nrst$subjectHits]  
    dist5p <- ifelse(abs(dist.s)<abs(dist.e), dist.s, dist.e) 
    
    ## get the lowest dist to either annot boundary from 3p side of the query
    dist.s <- end(query)[res.nrst$queryHits] - 
      start(subject)[res.nrst$subjectHits]    
    dist.e <- end(query)[res.nrst$queryHits] - 
      end(subject)[res.nrst$subjectHits]  
    dist3p <- ifelse(abs(dist.s)<abs(dist.e), dist.s, dist.e)     
  } else {
    #### no need to do calcs to start & end of subject since this clause assumes
    #### you have taken 5' or 3' of the subject!    
    ## get the lowest dist to annot boundary from 5p side of the query
    dist5p <- start(query)[res.nrst$queryHits] - 
      start(subject)[res.nrst$subjectHits]
    
    ## get the lowest dist to annot boundary from 3p side of the query
    dist3p <- end(query)[res.nrst$queryHits] - 
      start(subject)[res.nrst$subjectHits]    
  }
  
  ## get the lowest distance from the lowest 5p or 3p of the query!
  dist.lowest <- ifelse(abs(dist5p)<abs(dist3p), dist5p, dist3p) 
  
  ## fix signs to match biological upstream or downstream relative to query or subject! ##
  if(relativeTo=='query') {
    bore <- as.character(strand(query))[res.nrst$query]=="+"
    dist.lowest2 <- ifelse(bore, -dist.lowest, dist.lowest) 
  } else {
    bore <- as.character(strand(subject))[res.nrst$subjectHits]=="-"
    dist.lowest2 <- ifelse(bore, -dist.lowest, dist.lowest) 
  }
  rm(bore)
  
  res.nrst$dist <- dist.lowest2
  
  ## fix cases where two nested features were returned by choosing 
  ## the lowest absolute distances for both features.
  mins <- with(res.nrst, tapply(abs(dist), queryHits, min))
  res.nrst$lowest <- with(res.nrst, abs(dist)==mins[as.character(queryHits)])
  res.nrst <- droplevels(subset(res.nrst, lowest))
  res.nrst$lowest <- NULL
  
  res.nrst
}

#' Generate a window size label.
#'
#' Function to generate aesthetically pleasing window size label given an integer. This is one of the helper function used in \code{\link{getFeatureCounts}} & \code{\link{getFeatureCountsBig}}. 
#'
#' @param x vector of integers to generate the labels for. 
#'
#' @return a character vector of length(x) which has x normalized and suffixed by bp, Kb, Mb, or Gb depending on respective interval sizes.
#'
#' @seealso \code{\link{getFeatureCounts}}, \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' getWindowLabel(c(0,1e7,1e3,1e6,2e9))
getWindowLabel <- function(x) {
  ind <- cut(abs(x), c(0, 1e3, 1e6, 1e9, 1e12), 
             include.lowest = TRUE, right = FALSE, labels = FALSE)
  paste(x/c(1, 1e3, 1e6, 1e9, 1e12)[ind],
        c("bp", "Kb", "Mb", "Gb")[ind],sep="")
}

#' Get counts of annotation within a defined window around each query range/position. 
#'
#' Given a query object and window size(s), the function finds all the rows in subject which are <= window size/2 distance away. If weights are assigned to each positions in the subject, then tallied counts are multiplied accordingly. For large annotations, use \code{\link{getFeatureCountsBig}}.
#'
#' @param sites.rd RangedData/GRanges object to be used as the query.
#' @param features.rd RangedData/GRanges object to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves as a prefix to windows sizes!
#' @param chromSizes named vector of chromosome/space sizes to be used for testing if a position is off the mappable region. DEPRECATED and will be removed in future release.
#' @param widths a named/numeric vector of window sizes to be used for casting a net around each position. Default: \code{c(1000,10000,1000000)}.
#' @param weightsColname if defined, weigh each row from features.rd when tallying up the counts.
#' @param doInChunks break up sites.rd into small pieces of chunkSize to perform the calculations. Default is FALSE. Useful if you are expecting to find great deal of overlap between sites.rd and features.rd.
#' @param chunkSize number of rows to use per chunk of sites.rd. Default to 10000. Only used if doInChunks=TRUE.
#' @param parallel use parallel backend to perform calculation with \code{\link[foreach]{foreach}}. Defaults to FALSE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link[foreach]{registerDoSEQ}}.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd. There will be a column for each width defined in widths parameter. If widths was a named vector i.e. c("100bp"=100,"1K"=1000), then the colname parameter will be pasted together with width name else default name will be generated by the function.
#'
#' @note 
#' \itemize{
#'   \item If the input sites.rd parameter is GRanges object, then it is converted to RangedData and then converted back to GRanges at the end since \code{\link[IRanges]{findOverlaps}} function operates much faster on RangedData objects. 
#'   \item Try not to use this function for >50 spaces unless you have tons fo memory. 
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered before running the function. One can use any of the following libraries compatible with \code{\link[foreach]{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#' }
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}, \code{\link{getFeatureCountsBig}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene")
#' \dontrun{
#' geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene",doInChunks=TRUE, chunkSize=200)
#' geneCounts
#' ## Parallel version of getFeatureCounts
#' # geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene", parallel=TRUE)
#' # geneCounts
#' }
getFeatureCounts <- function(sites.rd, features.rd, 
                             colnam=NULL, chromSizes=NULL, 
                             widths=c(1000,10000,1000000), weightsColname=NULL, 
                             doInChunks=FALSE, chunkSize=10000, 
                             parallel=FALSE) {
  
  .checkArgsSetDefaults()
  
  if(!is.null(chromSizes)) {
    ## for historical version...display a warning if this parameter is passed!
    warning("decrepit option: chromSizes parameter is no longer required and will be ignored!")
  }    
  
  if(doInChunks & chunkSize < length(sites.rd)) {
  	rm("query","subject")
  	
    # no need to execute all this if chunkSize is bigger than data size!!!           
    total <- length(sites.rd)        
    starts <- seq(1,total,by=chunkSize)
    stops <- unique(c(seq(chunkSize ,total, by=chunkSize), total))
    stopifnot(length(starts)==length(stops))  
    message("Breaking up sites.rd into chunks of ",chunkSize)
    res <- GRanges()
    for(x in 1:length(starts)) {                
      res <- c(res,
               as(getFeatureCounts(sites.rd[starts[x]:stops[x],], 
                                   features.rd, colnam=colnam, widths=widths, 
                                   weightsColname=weightsColname, 
                                   parallel=parallel), "GRanges"))
    }

    if(RangedDataFlag) {
      res <- as(res,"RangedData")
    }

    return(res)
  } else {
    weighted <- ifelse(is.null(weightsColname),FALSE,TRUE)
    if(weighted) {
      mcols(subject)$weights <- mcols(features.rd)[,weightsColname]
    }
    rm(features.rd)
    
    # only get labels if not supplied
    if(is.null(names(widths))) {
      names(widths) <- getWindowLabel(widths)
    }            
    
    ## chunkize the objects for parallel processing ##
    chunks <- if(parallel) { 
      makeChunks(query, subject)
    } else {
      list(list("query"=query, "subject"=subject))
    }
        
    colnam <- cleanColname(colnam)
    
    ## perform overlap analysis in parallel by windows ##  
    res <- foreach(x=iter(chunks), .inorder=FALSE, .combine=rbind,
                   .export=c("widths","weighted","colnam"),
                   .packages=c("GenomicRanges","plyr")) %dopar% {                        
                     
                     counts <- sapply(widths, function(y) {
                       res.x <- findOverlaps(x$query, x$subject, 
                                             select='all', maxgap=(y/2), 
                                             ignore.strand=TRUE)
                       
                       if (weighted) {
                         res.x <- as.data.frame(res.x)
                         bore <- mcols(x$subject)[,"weights"][res.x$subjectHits]
                         res.x$weights <- bore
                         tapply(res.x$weights, res.x$queryHits, sum)            
                       } else {                                 
                         countQueryHits(res.x)
                       } 
                     })
                     counts <- as.data.frame(counts)
                     names(counts) <- paste(colnam,names(counts),sep=".")
                     counts$qID <- mcols(x$query)$tempyID
                     counts
                   }
    
    rm(chunks)
    
    ## merge results to the query object and return it ##
    .mergeAndReturn()
    
    sites.rd
  }
}

#' Get counts of annotation within a defined window around each query range/position for large annotation objects spanning greater than 1 billion rows.
#'
#' Given a query object and window size(s), the function finds all the rows in subject which are <= window size/2 distance away. Note that here counting is done using midpoint of the ranges in query instead of start-stop boundaries. The counts will differ slightly when compared to \code{\link{getFeatureCounts}}.
#'
#' @param sites.rd RangedData/GRanges object to be used as the query.
#' @param features.rd RangedData/GRanges object to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves as a prefix to windows sizes!
#' @param widths a named/numeric vector of window sizes to be used for casting a net around each position. Default: \code{c(1000,10000,1000000)}
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd. There will be a column for each width defined in widths parameter. If widths was a named vector i.e. c("100bp"=100,"1K"=1000), then the colname parameter will be pasted together with width name else default name will be generated by the function.
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}, \code{\link{getFeatureCounts}}.
#'
#' @export
getFeatureCountsBig <- function(sites.rd, features.rd, 
                                colnam=NULL, widths=c(1000,10000,1000000)) {
  
  .checkArgsSetDefaults()
  rm(features.rd)
   
  ranges(query) <- mid(ranges(query))
  query <- split(query, seqnames(query))
  subject <- split(subject,seqnames(subject)) 
  
  # only get labels if not supplied
  if(is.null(names(widths))) {
    names(widths) <- getWindowLabel(widths)
  }
    
  colnam <- cleanColname(colnam)
  
  ## get counts of midpoints using findInterval and add columns back to query object
  res <- lapply(names(widths), function(windowName) {
    message(".")
    columnName <- paste(colnam,names(widths[windowName]),sep=".")
        
    res.i <- lapply(ok.chrs, function(x) {    
      counts <- abs(findInterval(start(query[[x]]) - widths[windowName]/2 , 
                                 sort(start(subject[[x]]))) - 
                      findInterval(start(query[[x]]) + widths[windowName]/2, 
                                   sort(end(subject[[x]]))))
      res.x <- data.frame(qID=query[[x]]$tempyID)
      res.x[,columnName] <- counts
      res.x
    })
    
    do.call(rbind,res.i)    
  })
  
  res <- do.call(cbind,res)
  
  ## merge results to the query object and return it ##
  .mergeAndReturn()  
  
  sites.rd
}

#' Find overlapping positions/ranges that match between the query and subject. 
#'
#' When used in genomic context, the function annotates genomic positions of interest with information like if they were in a gene or cpg island or whatever annotation that was supplied in the subject.
#'
#' @param sites.rd RangedData/GRanges object to be used as the query.
#' @param features.rd RangedData/GRanges object to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves a core!
#' @param asBool Flag indicating whether to return results as TRUE/FALSE or the property of an overlapping feature..namely feature name and orientation if available. Defaults to FALSE.
#' @param feature.colnam column name from features.rd to be used for retrieving the feature name. By default this is NULL assuming that features.rd has a column that includes the word 'name' somewhere in it. Not required if asBool=TRUE.
#' @param parallel use parallel backend to perform calculation with \code{\link[foreach]{foreach}}. Defaults to FALSE. Not applicable when asBool=T. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link[foreach]{registerDoSEQ}}.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd.
#'
#' @note 
#' \itemize{
#'   \item If the input sites.rd parameter is GRanges object, then it is converted to RangedData and then converted back to GRanges at the end since \code{\link[IRanges]{findOverlaps}} function operates much faster on RangedData objects. 
#'   \item Try not to use this function for >50 spaces unless you have tons fo memory. 
#'   \item If parallel=TRUE, then be sure to have a parallel backend registered before running the function. One can use any of the following libraries compatible with \code{\link[foreach]{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doMC); registerDoMC(2)
#' }
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}}, \code{\link{getNearestFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene")
#' InGenes
#' \dontrun{
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene",asBool=TRUE)
#' InGenes
#' ## Parallel version of getSitesInFeature
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene",asBool=TRUE,parallel=TRUE)
#' InGenes
#' }
getSitesInFeature <- function(sites.rd, features.rd, colnam=NULL, 
                              asBool=FALSE, feature.colnam=NULL, 
                              parallel=FALSE) {    
  
  .checkArgsSetDefaults()
  
  ## chunkize the objects for parallel processing ##
  mcols(subject)$featureName <- mcols(features.rd)[,feature.colnam]
  rm(features.rd)
   
  chunks <- if(parallel) { 
    makeChunks(query, subject)
  } else {
    list(list("query"=query, "subject"=subject))
  }
    
  colnam <- cleanColname(colnam)
  
  ## perform overlap analysis in parallel by windows ##  
  res <- foreach(x=iter(chunks), .inorder=FALSE, .combine=rbind,
                 .export=c("colnam","asBool"),
                 .packages=c("GenomicRanges","plyr")) %dopar% {                        
                   
                   if (asBool) {
                     strand(x$subject) <- "*"
                     bore <- overlapsAny(x$query, x$subject, ignore.strand=TRUE)
                     res.x <- data.frame(qID=mcols(x$query)$tempyID, 
                                         featureName=bore)                   
                   } else {                                   
                     res.x <- as.data.frame(findOverlaps(x$query, x$subject, 
                                                         select='all', 
                                                         ignore.strand=TRUE))
                     res.x$qID <- mcols(x$query)$tempyID[res.x$queryHits]
                     
                     ## collapse rows where query returned two hits with the same featureNames 
                     ## due to alternative splicing or something else.    
                     res.x$featureName <- 
                       mcols(x$subject)[,"featureName"][res.x$subjectHits]
                     res.x$strand <- 
                       as.character(strand(x$subject))[res.x$subjectHits]					
                     
                     # isolate non-singletons to save time & memory! #
                     res.x <- merge(res.x, count(res.x,"queryHits"))
                     besties <- droplevels(subset(res.x,freq==1))
                     res.x <- droplevels(subset(res.x,freq>1))
                     
                     # use tapply instead of ddply() or by() because it's a lot faster on larger datasets #
                     bore <- with(res.x, 
                                  sapply(tapply(featureName,queryHits,unique), 
                                         paste, collapse=","))
                     res.x$featureName <- 
                       as.character(bore[as.character(res.x$queryHits)])
                     
                     bore <- with(res.x, sapply(tapply(strand,queryHits,unique), 
                                                paste, collapse=","))
                     res.x$strand <- 
                       as.character(bore[as.character(res.x$queryHits)])
                     
                     res.x <- unique(res.x[,c("queryHits","qID",
                                              "featureName","strand")])
                     
                     # put singletons & curated non-singletons back together! #
                     res.x <- rbind(besties[,names(res.x)], res.x)
                     rm(besties)
                     
                     res.x <- arrange(res.x, qID)
                     
                     names(res.x)[grepl("strand",names(res.x))] <- 
                       paste0(colnam,"Ort")
                   }
                   
                   ## change column names for swift merging by .mergeAndReturn() ##
                   names(res.x)[grepl("featureName",names(res.x))] <- colnam
                   
                   res.x
                 }
  
  rm(chunks)
  
  ## merge results to the query object and return it ##
  .mergeAndReturn()
  
  ## for legacy code support change sites.rd not in feature.rd to FALSE instead of NA ##
  if(is(sites.rd,"RangedData")) {
    sites.rd[[colnam]][is.na(sites.rd[[colnam]])] <- FALSE
  } else {
    mcols(sites.rd)[,colnam][is.na(mcols(sites.rd)[,colnam])] <- FALSE
  }
  
  sites.rd
}

#' Annotate a RangedData/GRanges object using one of annotation functions. 
#'
#' This is a wrapper function which calls one of following functions depending on annotType parameter: \code{\link{getFeatureCounts}}, \code{\link{getFeatureCountsBig}}, \code{\link{getNearestFeature}}, \code{\link{get2NearestFeature}}, code{\link{getSitesInFeature}} 
#'
#' @param annotType one of following: within, nearest, twoNearest, counts, countsBig.
#' @param ... Additional parameters to be passed to the respective annotation function.
#' @param postProcessFun function to call on the resulting object for any post processing steps.
#' @param postProcessFunArgs additional arguments for postProcessFun as a list.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd.
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}}, \code{\link{getFeatureCountsBig}}, \code{\link{getNearestFeature}}, \code{\link{get2NearestFeature}}, code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' doAnnotation(annotType="within",alldata.rd,genes.rd,"InGene",asBool=TRUE)
#' \dontrun{
#' doAnnotation(annotType="counts",alldata.rd,genes.rd,"NumOfGene")
#' doAnnotation(annotType="nearest",alldata.rd,genes.rd,"NearestGene")
#' doAnnotation(annotType="countsBig",alldata.rd,genes.rd,"ChipSeqCounts")
#' geneCheck <- function(x,wanted) { x$isWantedGene <- x$InGene %in% wanted; return(x) }
#' doAnnotation(annotType="within",alldata.rd,genes.rd,"InGene", postProcessFun=geneCheck, postProcessFunArgs=list("wanted"=c("FOXJ3","SEPT9","RPTOR")) )
#' }
doAnnotation <- function(annotType=NULL, ..., postProcessFun=NULL, 
                         postProcessFunArgs=list()) {    
  if(is.null(annotType)) {
    stop("Please define the annotType parameter to identify which type of annotation to perform: within, nearest, counts")
  }
  
  res <- switch(match.arg(annotType, c("within", "nearest", "twoNearest", 
                                       "counts", "countsBig")),
                within = getSitesInFeature(...),
                nearest = getNearestFeature(...),
                twoNearest = get2NearestFeature(...),
                counts = getFeatureCounts(...),
                countsBig = getFeatureCountsBig(...),     
                stop("Invalid annoType parameter")
  )
  
  if(!is.null(postProcessFun)) {
    res <- do.call(postProcessFun, append(postProcessFunArgs, 
                                          list(res), after=0))
  }
  
  res
}

#' Check args and set defaults.
#'
#' This function checks all the arguments passed to an annotation function and set default values for later use. Evaluation of this function happens in the parent function.
#'
.checkArgsSetDefaults <- function() {
  
  checks <- expression(
    RangedDataFlag <- FALSE,
    
    if(is(sites.rd,"RangedData")) {
      RangedDataFlag <- TRUE
      sites.rd <- as(sites.rd,"GRanges")
    },
    
    if(is(features.rd,"RangedData")) {
      features.rd <- as(features.rd,"GRanges")
    },
    
    if(!identical(class(sites.rd),class(features.rd))) {
      stop("sites.rd & features.rd are of different classes. 
           Please make them the same class: GRanges or RangedData")			
    },
    
    stopifnot(length(sites.rd)>0),
    stopifnot(length(features.rd)>0),
    
    if("parallel" %in% names(formals())) { 
      if(!parallel) { 
        registerDoSEQ() 
      }
    },
    
    if(exists("colnam")) {
      if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
      }
    },
    
    if (!any(unique(as.character(seqnames(sites.rd))) %in% 
               unique(as.character(seqnames(features.rd))))) {
      stop("There are no seqnames/spaces/chromosomes that are shared between the 
           query (sites.rd) and subject (features.rd)")
    },
    
    ## get feature names column for adding feature name to sites.rd ##
    ## dont throw an error if dists.only flag is TRUE from getNearestFeature ##
    if (exists("feature.colnam")) {
      if(is.null(feature.colnam)) {
        answer <- try(getRelevantCol(colnames(mcols(features.rd)),
                                     c("name","featureName"),
                                     "featureName",multiple.ok=TRUE), 
                      silent=TRUE)
        feature.colnam <- colnames(mcols(features.rd))[answer][1]
      }
      
      if(exists("dists.only")) {
        if(!dists.only & is.na(feature.colnam)) {
          stop("No featureName based column found.")
        }
      } else {
        if(is.na(feature.colnam)) { 
          stop("No featureName based column found.")
        }
      }			
    },
    
    ## check strand column of the feature ##
    if(all(strand(features.rd)=="*")) {
      message("Setting strand to '+' for features.rd parameter")
      strand(features.rd) <- "+"
    },
    
    ## check strand column of the feature ##
    if(all(strand(sites.rd)=="*")) {
      message("Setting strand to '+' for sites.rd parameter")
      strand(sites.rd) <- "+"
    },
    
    ## use only chromosomes that are present in both sites.rd and features.rd ##
    ok.chrs <- intersect(as.character(seqnames(sites.rd)),
                         as.character(seqnames(features.rd))),
    features.rd <- keepSeqlevels(features.rd, ok.chrs),
    good.rows <- as.character(seqnames(sites.rd)) %in% ok.chrs,    
    
    ## extract required objects to streamline downstream code/steps ##
    ## tag each row with tempyID for merging with original object ##
    ## this is crucial since objects are divided into chunks which resets the index from 1...n ##
    ## tempyID would preserve the original order for parallel processing ##
    query <- NULL,
    query <- sites.rd,
    mcols(query) <- NULL,
    mcols(query)$tempyID <- 1:length(query),
    mcols(sites.rd)$tempyID <- mcols(query)$tempyID,
    
    subject <- NULL,
    subject <- features.rd,
    mcols(subject) <- NULL,
    mcols(subject)$tempyID <- 1:length(subject)
  )
  
  eval.parent(checks)
}

#' Merge results back to the query object and perform additional post processing  steps.
#'
#' This function merges all the calculation results back to the query object. Additionally, if any flags were set, the function does the necessary checks and processing to format the return object as required. Evaluation of this function happens in the parent function.
#'
.mergeAndReturn <- function() {
  toDo <- expression(
    
    ## make sure res object has the required fields ##
    stopifnot("qID" %in% names(res)),
    
    ## make sure we only have one resulting row per query ##
    stopifnot(!any(table(res$qID)>1)),
    
    res <- DataFrame(res),
    
    ## setup new columns to be added using NA and add the proper class ## 
    newCols <- grep(colnam,names(res),value=TRUE),
    mcols(sites.rd)[newCols] <- NA,
    for(newCol in newCols) {
      mcols(sites.rd)[[newCol]] <- as(mcols(sites.rd)[[newCol]], 
                                      class(res[[newCol]]))
    },
    
    ## merge back results in the same order as rows in query/sites.rd ##
    rows <- match(mcols(sites.rd)$tempyID[good.rows], res$qID),
    mcols(sites.rd)[good.rows,][!is.na(rows),newCols] <- 
      res[rows[!is.na(rows)], newCols],
    
    ## clear up any temp columns ##
    mcols(sites.rd)$tempyID <- NULL,
    
    if(RangedDataFlag) {
      sites.rd <- as(sites.rd,"RangedData")
    }
  )
  
  eval.parent(toDo)
}
