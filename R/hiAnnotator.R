#' Annotating RangedData or GRanges objects with hiAnnotator.
#'
#' hiAnnotator contains set of functions which allow users to annotate a RangedData or GRanges object with custom set of annotations. The basic philosophy of this package is to take two RangedData or GRanges objects (query & subject) with common set of space/seqnames (i.e. chromosomes) and return associated annotation per space and rows from the query matching space and rows from the subject (i.e. genes or cpg islands).
# The package comes with three types of annotation functions which calculates if a position from query is: within a feature, near a feature, or count features in defined window sizes. Moreover, one can utilize parallel backend for each annotation function to utilize multiple processors. In addition, the package is equipped with a wrapper function, which finds appropriate columns needed to make a RangedData or GRanges object from a common dataframe..
#'
#' @import IRanges foreach iterators doBy RMySQL rtracklayer dataframe GenomicRanges BSgenome
#' @docType package
#' @name hiAnnotator
NULL

#' Sample HIV Integration Sites data
#' 
#' A sample dataset containing collection of unique HIV integration sites in the human genome mapped to UCSC freeze hg18.
#' 
#' \itemize{
#'   \item id. A unique id for each row.
#'   \item Position. The genomic coordinate of the integration site.
#'   \item Chr. The chromosome of the integration site. 
#'   \item Ort. The orientation or strand of the integration site. 
#'   \item sampleName. Patient from which the gDNA sample was taken (artificially generated).
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 2408 rows and 5 variables
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
#' @seealso \code{\link{getUCSCtable}}, \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' session <- makeUCSCsession()
#' genome(session)
#' session <- makeUCSCsession("mm8")
#' genome(session)
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
#' refflat <- getUCSCtable("refFlat","RefSeq Genes") ## same as session <- makeUCSCsession(); 
#' refflat <- getUCSCtable("refFlat", "RefSeq Genes", bsession=session, freeze="hg18")
#' head(refflat)
getUCSCtable <- function(tableName, trackName, bsession=NULL, freeze="hg18", ...) {
    if(is.null(bsession)) { 
        bsession <- makeUCSCsession(freeze)
    }
    if(!tableName %in% tableNames(ucscTableQuery(bsession,track=trackName))) {
        stop(paste("The provided table name:",tableName,"doesn't exists in track",trackName,"on UCSC for",freeze,"genome"))
    }
    
    ## using getTable() instead of track() due to "No supported output types" error for certain annotation types.    
    getTable(ucscTableQuery(bsession,track=trackName,table=tableName,...))
}

#' Find the column index of interest given the potential choices.
#'
#' The function finds relevant column(s) of interest from a vector of column names derived from a dataframe or a RangedData object. If no usable column is found, the function spits out a relevant error or returns the index of the usable column(s). This is an assistant function called by functions listed in the see also section.
#'
#' @param col.names column names from a dataframe or a RangedData object
#' @param col.options potential column names or partial names that may exist in col.names
#' @param col.type type of column information the function is searching for, used in construction of error messages. 
#' @param multiple.ok if multiple matches are found then return indices, else spit an error out.
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
getRelevantCol <- function(col.names, col.options, col.type, multiple.ok=FALSE) {
    answer <- unique(as.numeric(unlist(sapply(col.options,function(x) grep(x,col.names,ignore.case=TRUE)))))
    if(length(answer)>1) {
        if(!multiple.ok) {
            stop(paste("More than one",col.type,"based column found:",paste(col.names[answer],sep="",collapse=", ")))
        } else {
            answer
        }
    } else if (length(answer)==0) {
        stop(paste("No",col.type,"based column found."))
    } else {
        answer
    }
}

#' Make a sorted RangedData object from a dataframe. 
#'
#' The function converts a dataframe into a RangedData object without too much hassle of renaming column names. The function finds column names that sound like space, chromosome, start, stop, position, etc and puts them in respective slots to facilitate the conversion of a dataframe to RangedData object. If more than one column that sounds like start, stop, or position is present, the function will use the first match as the representative. It is recommended to run this function before utilizing any other annotation functions since it will sort the object by chromosome and position for copying annotations back to their respective rows accurately.
#'
#' @param x dataframe to be converted into a RangedData object
#' @param positionsOnly boolean flag indicating to return only position based data or everything from the dataframe. Defaults to FALSE.
#' @param soloStart flag denoting whether only one position based column is available. In other words, only starts are present and no stops. Default=FALSE.
#' @param chromCol use the defined column name for space/chromosome based data from the dataframe. Defaults to NULL.
#' @param strandCol use the defined column name for strand or orientation from the dataframe. Defaults to NULL.
#' @param startCol use the defined column name for start coordinate from the dataframe. Defaults to NULL.
#' @param stopCol use the defined column name for stop coordinate from the dataframe. Defaults to NULL and not required if soloStart=TRUE.
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
#' makeRangedData(sites,soloStart=TRUE)
#' makeRangedData(sites,soloStart=TRUE,positionsOnly=TRUE)
#' makeRangedData(sites)
#'
#' data(genes)
#' head(genes)
#'
#' makeRangedData(genes,soloStart=TRUE)
#' makeRangedData(genes)
makeRangedData <- function(x, positionsOnly=FALSE, soloStart=FALSE, chromCol=NULL, strandCol=NULL, startCol=NULL, stopCol=NULL) {
    ## set column names for space and strand if not provided ##
    if(is.null(chromCol)) {
        colIndex <- getRelevantCol(names(x),c("chr","chromosome","tname","space","chrom","contig","seqnames"),"space")
        names(x)[colIndex] <- "space"
    } else {
        names(x)[names(x)==chromCol] <- "space"
    }
    
    if(is.null(strandCol)) {
        colIndex <- getRelevantCol(names(x),c("ort","orientation","strand"),"strand")
        names(x)[colIndex] <- "strand"
    } else {
        names(x)[names(x)==strandCol] <- "strand"
    }    
        
    if(is.null(startCol)) {
        startCol <- getRelevantCol(names(x),c("position","intsite","txstart","start","chromstart"),"start",multiple.ok=TRUE)
        names(x)[startCol[1]] <- "start" 
    } else {
        names(x)[names(x)==startCol] <- "start"
    }
    
    ## only do stop if soloStart=F ##
    if(!as.logical(soloStart) & is.null(stopCol)) {
        stopCol <- getRelevantCol(names(x),c("txend","end","stop","chromend"),"end",multiple.ok=TRUE)
        names(x)[stopCol[1]] <- "end" 
    } else {
        names(x)[names(x)==stopCol] <- "end"
    }
    
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
        x <- orderBy(~space+mid,x)
        x$mid <- NULL
    } else {  
        x <- orderBy(~space+start,x)
        x$end <- x$start
    }
    
    if(as.logical(positionsOnly)) {
        sites.rd <- as(x[,c("space","start","end","strand")],"RangedData")
    } else {
        sites.rd <- as(x,"RangedData")
    }
    
    sites.rd    
}

#' Make a sorted GRanges object from a dataframe. 
#'
#' The function converts a dataframe into a GRanges object without too much hassle of renaming column names. The function finds column names that sound like seqname, chromosome, start, stop, position, etc and puts them in respective slots to facilitate the conversion of a dataframe to GRanges object. If more than one column that sounds like start, stop, or position is present, the function will use the first match as the representative. It is recommended to run this function before utilizing any other annotation functions since it will sort the object by chromosome and position for copying annotations back to their respective rows accurately. This function wraps around \code{\link{makeRangedData}}.
#'
#' @param x dataframe to be converted into a GRanges object
#' @param freeze UCSC genome version of the data in x. Default is NULL. This parameter is generally used to populate \code{\link{seqinfo}} slot of GRanges objects.
#' @param ... parameters for \code{\link{makeRangedData}}.
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
#' makeGRanges(genes,soloStart=TRUE)
#' makeGRanges(genes, freeze="hg18", soloStart=TRUE)
#' makeGRanges(genes)
#' makeGRanges(genes, freeze="hg18")
makeGRanges <- function(x, freeze=NULL, ...) {    
	sites.gr <- as(makeRangedData(x, ...), "GRanges")	
	
	if(!is.null(freeze)) {
		genomeLib <- grep(freeze,installed.genomes(),value=TRUE)
		if(length(genomeLib)!=0) {		
			bsGenomeObject <- strsplit(genomeLib,"\\.")[[1]][2]
			chrom.info <- seqlengths(do.call(`:::`,list(genomeLib,bsGenomeObject)))
		} else {
			## get the chromInfo file from UCSC
			z <- gzcon(url(paste0("http://hgdownload.cse.ucsc.edu/goldenPath/", freeze, "/database/chromInfo.txt.gz")))
			zlines <- try(readLines(z))
			close(z)
			if (class(zlines)=="try-error") stop("Could not get thru to UCSC server - try later!")
			raw.data <- textConnection(zlines)
			chrom.info <- read.delim(raw.data,header=F,stringsAsFactors=F)[,1:2]
			chrom.info <- structure(chrom.info$V2, names=chrom.info$V1)
			close(raw.data)
		}
	
		# append to sites.gr as Seqinfo slot
		chrom.info <- chrom.info[names(chrom.info) %in% levels(seqnames(sites.gr))]
		chrom.info <- chrom.info[order(suppressWarnings(as.numeric(gsub("chr","",names(chrom.info)))))] ## reorder things while at it!
		seqinfo(sites.gr,new2old=order(suppressWarnings(as.numeric(gsub("chr","",levels(seqnames(sites.gr))))))) <- Seqinfo(names(chrom.info), chrom.info, rep(FALSE,length(chrom.info)), freeze)
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
#' @param side boundary of annotation to use to calculate the nearest distance. Options are '5p','3p', or the default 'either'.
#' @param feature.colnam column name from features.rd to be used for retrieving the nearest feature name. By default this is NULL assuming that features.rd has a column that includes the word 'name' somewhere in it.
#' @param strand.colnam column name from features.rd to be used for retrieving the nearest feature's orientation. By default this is NULL assuming that features.rd has a column that includes the word 'strand' somewhere in it. If it doesn't the function will assume the supplied annotation is in '+' orientation (5' -> 3'). The same applies to strand column in sites.rd.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to FALSE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd.
#'
#' @note Try not to use this function for >50 spaces unless you have tons fo memory. If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}, \code{\link{get2NearestFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' head(sites)
#' alldata.rd <- makeRangedData(sites,soloStart=TRUE)
#' alldata.rd
#'
#' data(genes)
#' head(genes)
#' genes.rd <- makeRangedData(genes)
#' genes.rd
#'
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene")
#' nearestGenes
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",side="5p")
#' nearestGenes
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene",side="3p")
#' nearestGenes
#'
#' # Parallel version of getNearestFeature
#' nearestGenes <- getNearestFeature(alldata.rd,genes.rd,"NearestGene", parallel=TRUE)
#' nearestGenes
getNearestFeature <- function(sites.rd, features.rd, colnam=NULL, side="either", feature.colnam=NULL, strand.colnam=NULL, parallel=FALSE) {
    stopifnot(nrow(sites.rd)>0)
    stopifnot(nrow(features.rd)>0)
    
    grangesFlag <- FALSE
    if(class(sites.rd)=="GRanges" & class(features.rd)=="GRanges") {
    	grangesFlag <- TRUE
    	sites.rd <- as(sites.rd,"RangedData")
    	features.rd <- as(features.rd,"RangedData")
    }
    
    if (!any(names(sites.rd) %in% names(features.rd))) {
        stop("There are no spaces/chromosomes that are shared between the query (sites.rd) and subject (features.rd)")
    }
    
    if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
    }
    
    ## use only chromosomes that are present in both sites.rd and features.rd ##
    ok.chrs <- intersect(space(sites.rd),space(features.rd))
    features.rd <- features.rd[names(features.rd) %in% ok.chrs]
    
    if(is.null(feature.colnam)) {
        featureName <- getRelevantCol(colnames(features.rd),c("name","featureName"),"featureName",multiple.ok=TRUE)
        feature.colnam <- colnames(features.rd)[featureName][1]
        message("Using column ", feature.colnam, " for feature.colnam parameter.")
    }
    
    ## do a check of strand column in sites.rd ##
    answer <- unique(as.numeric(unlist(sapply(c("strand","orientation","ort"), function(x) grep(x,colnames(sites.rd), ignore.case = TRUE)))))
    strandCol <- colnames(sites.rd)[answer][1]
    if(is.na(strandCol)) { 
        message("No orientation column found in sites.rd. Using '+' as default.")
        sites.rd$strand <- "+" 
    }
    
    ## do a check of strand column in features.rd ##
    answer <- unique(as.numeric(unlist(sapply(c("strand","orientation","ort"), function(x) grep(x,colnames(features.rd), ignore.case = TRUE)))))
    strandCol <- colnames(features.rd)[answer][1]
    if(is.null(strand.colnam) & is.na(strandCol)) { 
        message("No orientation column found in features.rd. Using '+' as default.")
        features.rd$strand <- "+"
    }
    if(!is.na(strandCol) & strandCol!="strand") { features.rd$strand <- unlist(values(features.rd))[strandCol] }
    
    ## convert any factor columns to character to avoid downstream issues with NAs and unlisting of CompressedCharacterList object ##
    factorCols <- sapply(colnames(features.rd),function(x) class(features.rd[[x]]))=="factor"
    if(any(factorCols)) {
        for (x in names(which(factorCols))) { features.rd[[x]] <- as.character(features.rd[[x]]) }        
    }
    
    ## extract required objects to streamline downstream code ##
    query <- ranges(sites.rd)

    subject <- ranges(features.rd) ## start with both sides of features...aka side='either'
    vals.s <- values(features.rd)
    
    if(side %in% c('5p','3p')) {
        ##get only 5 prime sides of features
        if (side=='5p')
            subject <- ranges(RangedData(IRanges(start=ifelse(features.rd$strand=="+", start(features.rd), end(features.rd)), width=1), 
            							 space=space(features.rd)))  
        
        ##get only 3 prime sides of features
        if (side=='3p')
            subject <- ranges(RangedData(IRanges(start=ifelse(features.rd$strand=="-", start(features.rd), end(features.rd)), width=1),
            							 space=space(features.rd)))                     
    }
    
    if(!parallel) { registerDoSEQ() }
    
    ## first get the nearest indices ##
    res <- foreach(x=iter(ok.chrs), .inorder=TRUE, .export=c("query","subject"), .packages="IRanges") %dopar% { 
		as.matrix(nearest(query[[x]], subject[[x]], select="all")) 
    }     
    names(res) <- ok.chrs
    stopifnot(identical(lapply(query,length)[ok.chrs],lapply(res,function(x) length(unique(x[,"queryHits"])) )[ok.chrs])) ## check for safety
    
    # check if >1 nearest matches found, get the indices of the feature with shortest distance to 5p/3p
    res <- getLowestDists(query,subject,subjectOrt=vals.s[,"strand"],ok.chrs=ok.chrs,res.nrst=res,side=side,cores.use=1)
    
    ## fix any cases where matrix got converted to integer due to only value found ##
    if(any(unlist(lapply(res,class))!="matrix")) {
        tofix <- names(which(unlist(lapply(res,class))!="matrix"))
        for (i in tofix) { res[[i]] <- t(res[[i]]) }
    }
    
    ## for the feature of shortest indices, get the names, and strand attributes
    featureName <- foreach(x=iter(ok.chrs), .inorder=TRUE, .export=c("vals.s","res","feature.colnam")) %dopar% {
		vals.s[[x]][res[[x]][,"subjectHits"],feature.colnam] 
    }     
    
    ort <- foreach(x=iter(ok.chrs), .inorder=TRUE, .export=c("vals.s","res","feature.colnam")) %dopar% { 
    	vals.s[[x]][res[[x]][,"subjectHits"],"strand"]
    }     
    
    names(featureName) <- names(ort) <- ok.chrs
    
    stopifnot(identical(lapply(res,nrow),lapply(featureName,length))) ## check for safety
    stopifnot(identical(lapply(res,nrow),lapply(ort,length))) ## check for safety
    
    ## fix cases where two equally nearest features were returned by concatenating feature names and Ort while returning one distance per query
    res.i <- foreach(x=iter(ok.chrs), .inorder=TRUE, .export=c("ort","res","featureName")) %dopar% {
        toprune <- as.data.frame(res[[x]])
        toprune$featureName <- featureName[[x]]
        toprune$ort <- ort[[x]]
        toprune <- unique(toprune[,-2])
        counts <- table(toprune$queryHits)
        ismulti <- toprune$queryHits %in% as.numeric(names(counts[counts>1]))
        if(any(ismulti)) {
            goods <- toprune[!ismulti,]
            toprune <- toprune[ismulti,]
            tempStore <- with(toprune,sapply(tapply(featureName,queryHits,unique),paste,collapse=","))
            toprune$featureName <- as.character(tempStore[as.character(toprune$queryHits)])            
            tempStore <- with(toprune,sapply(tapply(ort,queryHits,unique),paste,collapse=","))
            toprune$ort <- as.character(tempStore[as.character(toprune$queryHits)])            
            tempStore <- with(toprune,sapply(tapply(lowestDist,queryHits,abs),min)) ## if a site falls exactly between two genes pick one the abs(lowest)
            toprune$lowestDist <- as.numeric(tempStore[as.character(toprune$queryHits)])
            toprune <- rbind(goods,unique(toprune))
        }
        return(orderBy(~queryHits,toprune))
    }  
    names(res.i) <- ok.chrs
    
    ## add columns back to query object
    prefix <- ifelse(side=="either","",side)
    good.rows <- space(sites.rd) %in% ok.chrs
    colnam <- cleanColname(colnam)
    
    sites.rd[[paste(prefix,colnam,sep="")]][good.rows] <- unsplit(lapply(res.i,"[[","featureName"), space(sites.rd)[good.rows], drop = TRUE)
    sites.rd[[paste(prefix,colnam,"Ort",sep="")]][good.rows] <- unsplit(lapply(res.i,"[[","ort"), space(sites.rd)[good.rows], drop = TRUE)
    sites.rd[[paste(prefix,colnam,"Dist",sep="")]][good.rows] <- unsplit(lapply(res.i,"[[","lowestDist"), space(sites.rd)[good.rows], drop = TRUE)    
    
    if(grangesFlag) {
    	sites.rd <- as(sites.rd,"GRanges")
    }
    
    sites.rd
}

#' Get two nearest upstream and downstream annotation boundary for a position range. 
#'
#' Given a query object, the function retrieves the two nearest feature upstream and downstream along with their properties from a subject and then appends them as new columns within the query object. When used in genomic context, the function can be used to retrieve two nearest gene upstream and downstream of the genomic position of interest.
#'
#' @param sites.rd RangedData/GRanges object to be used as the query.
#' @param features.rd RangedData/GRanges object to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves a core!
#' @param side boundary of annotation to use to calculate the nearest distance. Options are '5p','3p', or the default 'either'.
#' @param feature.colnam column name from features.rd to be used for retrieving the nearest feature name. By default this is NULL assuming that features.rd has a column that includes the word 'name' somewhere in it.
#' @param strand.colnam column name from features.rd to be used for retrieving the nearest feature's orientation. By default this is NULL assuming that features.rd has a column that includes the word 'strand' somewhere in it. If it doesn't the function will assume the supplied annotation is in '+' orientation (5' -> 3'). The same applies to strand column in sites.rd.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd.
#'
#' @note For cases where a position is at the edge and there are no feature up/down stream since it would fall off the chromosome, the function simply returns the nearest feature. In addition, if there are multiple locations where a query falls into, the function arbitrarily chooses one to serve as the nearest feature, then reports 2 upstream & downstream feature. That may occasionally yield features which are the same upstream and dowstream, which is commonly encountered when studying spliced genes or phenomena related to it. 
#'
#' @seealso \code{\link{getNearestFeature}}, \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' head(sites)
#' alldata.rd <- makeRangedData(sites,soloStart=TRUE)
#' alldata.rd
#'
#' data(genes)
#' head(genes)
#' genes.rd <- makeRangedData(genes)
#' genes.rd
#'
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene")
#' nearestGenes
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene",side="5p")
#' nearestGenes
#' nearestGenes <- get2NearestFeature(alldata.rd,genes.rd,"NearestGene",side="3p")
#' nearestGenes
#'
get2NearestFeature <- function(sites.rd, features.rd, colnam=NULL, side="either", feature.colnam=NULL, strand.colnam=NULL) {
	stopifnot(nrow(sites.rd)>0)
    stopifnot(nrow(features.rd)>0)
    
    grangesFlag <- FALSE
    if(class(sites.rd)=="GRanges" & class(features.rd)=="GRanges") {
    	grangesFlag <- TRUE
    	sites.rd <- as(sites.rd,"RangedData")
    	features.rd <- as(features.rd,"RangedData")
    }
    
    if (!any(names(sites.rd) %in% names(features.rd))) {
        stop("There are no spaces/chromosomes that are shared between the query (sites.rd) and subject (features.rd)")
    }
    
    if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
    }
    
    ## use only chromosomes that are present in both sites.rd and features.rd ##
    ok.chrs <- intersect(space(sites.rd),space(features.rd))
	features.rd <- features.rd[names(features.rd) %in% ok.chrs]    

    if(is.null(feature.colnam)) {
        featureName <- getRelevantCol(colnames(features.rd),c("name","featureName"),"featureName",multiple.ok=TRUE)
        feature.colnam <- colnames(features.rd)[featureName][1]  
    }
    
    ## do a check of strand column in sites.rd ##
    answer <- unique(as.numeric(unlist(sapply(c("strand","orientation","ort"), function(x) grep(x,colnames(sites.rd), ignore.case = TRUE)))))
    strandCol <- colnames(sites.rd)[answer][1]
    if(is.na(strandCol)) { 
        message("No orientation column found in sites.rd. Using '+' as default.")
        sites.rd$strand <- "+" 
    }
    
    ## do a check of strand column in features.rd ##
    answer <- unique(as.numeric(unlist(sapply(c("strand","orientation","ort"), function(x) grep(x,colnames(features.rd), ignore.case = TRUE)))))
    strandCol <- colnames(features.rd)[answer][1]
    if(is.null(strand.colnam) & is.na(strandCol)) { 
        message("No orientation column found in features.rd. Using '+' as default.")
        features.rd$strand <- "+"
    }
    if(!is.na(strandCol) & strandCol!="strand") { features.rd$strand <- unlist(values(features.rd))[strandCol] }

    ## convert any factor columns to character to avoid downstream issues with NAs and unlisting of CompressedCharacterList object ##
    factorCols <- sapply(colnames(features.rd),function(x) class(features.rd[[x]]))=="factor"
    if(any(factorCols)) {
        for (x in names(which(factorCols))) { features.rd[[x]] <- as.character(features.rd[[x]]) }        
    }
    
    ## extract objects/lists of interest for dowstream code
    vals.q <- values(sites.rd) 
    vals.s <- values(features.rd)
    
    query <- ranges(sites.rd)

    subject <- ranges(features.rd) ## start with both sides of features...aka side='either'
    
    if(side %in% c('5p','3p')) {
        ##get only 5 prime sides of features
        if (side=='5p')
            subject <- ranges(RangedData(IRanges(start=ifelse(features.rd$strand=="+",start(features.rd),end(features.rd)),width=1),
            							 space=space(features.rd)))  
        
        ##get only 3 prime sides of features
        if (side=='3p')
            subject <- ranges(RangedData(IRanges(start=ifelse(features.rd$strand=="-",start(features.rd),end(features.rd)),width=1),
            							 space=space(features.rd)))                     
    }
    
    ## u = upstream, d = downstream
    ## thinking concept: u2.....u1.....intSite(+).....d1.....d2
    ## thinking concept: d2.....d1.....intSite(-).....u1.....u2
    ## searching concept: res.left2.....res.left1.....res....intSite....res.....res.right1.....res.right2
    res <- sapply(ok.chrs,function(x) nearest(query[[x]],subject[[x]]))
    res.left1 <- sapply(ok.chrs,function(x) res[[x]]-1)
    res.left2 <- sapply(ok.chrs,function(x) res[[x]]-2) 
    res.right1 <- sapply(ok.chrs,function(x) res[[x]]+1)
    res.right2 <- sapply(ok.chrs,function(x) res[[x]]+2) 
        
    dist.nrst <- getLowestDists(query,subject,vals.s[,"strand"],ok.chrs,res,side,1)
    
    ## perform upstream-downstream checks by testing distances
    u1 <- sapply(ok.chrs, function(x) ifelse(dist.nrst[[x]]<0, res[[x]], ifelse(vals.q[[x]][,"strand"]=="+",res.left1[[x]],res.right1[[x]])))   
    d1 <- sapply(ok.chrs, function(x) ifelse(dist.nrst[[x]]<0, ifelse(vals.q[[x]][,"strand"]=="+",res.right1[[x]],res.left1[[x]]), res[[x]]))
    u2 <- sapply(ok.chrs, function(x) ifelse(dist.nrst[[x]]<0, ifelse(vals.q[[x]][,"strand"]=="+",res.left1[[x]],res.right1[[x]]), ifelse(vals.q[[x]][,"strand"]=="+",res.left2[[x]],res.right2[[x]])))
    d2 <- sapply(ok.chrs, function(x) ifelse(dist.nrst[[x]]<0, ifelse(vals.q[[x]][,"strand"]=="+",res.right2[[x]],res.left2[[x]]), ifelse(vals.q[[x]][,"strand"]=="+",res.right1[[x]],res.left1[[x]])))

    prefix <- ifelse(side=="either","Either",side)
    good.rows <- space(sites.rd) %in% ok.chrs
    
    message("u = upstream, d = downstream")
    message("thinking concept: u2.....u1.....intSite(+).....d1.....d2")
    message("thinking concept: d2.....d1.....intSite(-).....u1.....u2")

    colnam <- cleanColname(colnam)

    ## add columns back to query object
    for (f in c("u1","u2","d1","d2")) {
        message(f)
        res.nrst <- get(f)
        res.nrst <- sapply(ok.chrs, function(x) ifelse(res.nrst[[x]]>length(subject[[x]]) | res.nrst[[x]]<1,res[[x]],res.nrst[[x]])) 
        features <- sapply(ok.chrs, function(x) vals.s[[x]][res.nrst[[x]],feature.colnam])
        ort <- sapply(ok.chrs, function(x) vals.s[[x]][res.nrst[[x]],"strand"]) 
        dists <- getLowestDists(query,subject,vals.s[,"strand"],ok.chrs,res.nrst,side,2)
        
        if (f == "u1") { coldef <- paste(prefix,colnam,"upStream1",sep=".") }
        if (f == "u2") { coldef <- paste(prefix,colnam,"upStream2",sep=".") }
        if (f == "d1") { coldef <- paste(prefix,colnam,"downStream1",sep=".") }
        if (f == "d2") { coldef <- paste(prefix,colnam,"downStream2",sep=".") }
        
        sites.rd[[coldef]][good.rows] <- unsplit(features, space(sites.rd)[good.rows], drop = TRUE)
        sites.rd[[paste(coldef,"Ort",sep=".")]][good.rows] <- unsplit(ort, space(sites.rd)[good.rows], drop = TRUE)
        sites.rd[[paste(coldef,"Dist",sep=".")]][good.rows] <- unsplit(dists, space(sites.rd)[good.rows], drop = TRUE)    
    }
    
    if(grangesFlag) {
    	sites.rd <- as(sites.rd,"GRanges")
    }
            
    sites.rd
}

#' Get the lowest distance from the 5' or 3' boundaries of query and subject.
#'
#' Given a query and subject with indicies from \code{\link{nearest}}, calculate the shortest distance to either boundaries of the query and subject. This is a helper function utilized in \code{\link{getNearestFeature}}, \code{\link{get2NearestFeature}}
#'
#' @param query Compressed RangedList object to be used as the query.
#' @param subject Compressed RangedList object to be used as the subject.
#' @param subjectOrt CompressedCharacterList of orientation/strand of features in subject broken up by ok.chrs...generated using values(subject)[,"strand"]
#' @param ok.chrs vector of spaces or chromosomes present in query/subject to do the calculations for.
#' @param res.nrst nearest indices as returned by \code{\link{nearest}}.
#' @param side boundary of subject/annotation to use to calculate the nearest distance. Options are '5p','3p', or the default 'either'.
#' @param cores.use number of cores to utilize for the calculation. Default is 1. If >1, then a parallel backend must be registered to perform calculation with \code{\link{foreach}}.
#'
#' @return a list of lowest distances named by ok.chrs
#'
#' @seealso \code{\link{getNearestFeature}}, \code{\link{get2NearestFeature}}.
#'
#' @export
getLowestDists <- function(query=NULL, subject=NULL, subjectOrt=NULL, ok.chrs=NULL, res.nrst=NULL, side="either", cores.use=1) {
    if(is.null(query) | is.null(subjectOrt) | is.null(subject) | is.null(ok.chrs) | is.null(res.nrst)) {
        stop("One of following is null: query, subjectOrt, subject, ok.chrs, res.nrst")
    }    
    
    if(any(!names(subject) %in% names(subjectOrt))) {
        stop("There are missing spaces/chromosomes in subjectOrt which are present in subject")
    }
    
    ## convert res.nrst to a matrix if its a list of integers/numeric resulting from using the default nearest() used in get2NearestFeature()
    ismatrixFlag <- TRUE
    if(all(unlist(lapply(res.nrst,class))!="matrix")) {
        res.nrst <- lapply(res.nrst,function(x) as.matrix(cbind(queryHits=1:length(x),subjectHits=x)) )
        ismatrixFlag <- FALSE
    }

    if(any(unlist(sapply(ok.chrs,function(x) res.nrst[[x]][,"subjectHits"]>length(subject[[x]]) | res.nrst[[x]][,"subjectHits"]<1)))) {
        stop("Subject indicies out of vector range in res.nrst ... this can introduce NAs and cause science to fail!")
    }
    
    if(cores.use==1) { registerDoSEQ() }
    
    if (side=="either") {
        ## get the lowest dist to either annot boundary from 5p side of the query
        dist.s <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { start(query[[x]][res.nrst[[x]][,"queryHits"]])-start(subject[[x]][res.nrst[[x]][,"subjectHits"]]) }     
        dist.e <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { start(query[[x]][res.nrst[[x]][,"queryHits"]])-end(subject[[x]][res.nrst[[x]][,"subjectHits"]]) }     
        names(dist.s) <- names(dist.e) <- ok.chrs

        dist5p <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("dist.s","dist.e")) %dopar% { ifelse(abs(dist.s[[x]])<abs(dist.e[[x]]),dist.s[[x]],dist.e[[x]]) }

        ## get the lowest dist to either annot boundary from 3p side of the query
        dist.s <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { end(query[[x]][res.nrst[[x]][,"queryHits"]])-start(subject[[x]][res.nrst[[x]][,"subjectHits"]]) }     
        dist.e <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { end(query[[x]][res.nrst[[x]][,"queryHits"]])-end(subject[[x]][res.nrst[[x]][,"subjectHits"]]) }     
        names(dist.s) <- names(dist.e) <- ok.chrs

        dist3p <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("dist.s","dist.e")) %dopar% { ifelse(abs(dist.s[[x]])<abs(dist.e[[x]]),dist.s[[x]],dist.e[[x]]) }               
    } else {
        ## get the lowest dist to annot boundary from 5p side of the query
        dist5p <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { start(query[[x]][res.nrst[[x]][,"queryHits"]])-start(subject[[x]][res.nrst[[x]][,"subjectHits"]]) }
        
        ## get the lowest dist to annot boundary from 3p side of the query
        dist3p <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { end(query[[x]][res.nrst[[x]][,"queryHits"]])-start(subject[[x]][res.nrst[[x]][,"subjectHits"]]) }
    }
    
    names(dist3p) <- names(dist5p) <- ok.chrs
        
    ## get the lowest distance from the lowest 5p or 3p of the query!
    dist.lowest <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("dist5p","dist3p")) %dopar% { ifelse(abs(dist5p[[x]])<abs(dist3p[[x]]),dist5p[[x]],dist3p[[x]]) }
    names(dist.lowest) <- ok.chrs
    
    stopifnot(identical(unlist(lapply(res.nrst[ok.chrs],nrow)),unlist(lapply(dist.lowest,length)))) ## check for safety
    
    ## fix signs to match biological upstream or downstream
    dist.lowest2 <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("subjectOrt","res.nrst","dist.lowest")) %dopar% { ifelse(subjectOrt[[x]][res.nrst[[x]][,"subjectHits"]]=="-",-dist.lowest[[x]],dist.lowest[[x]]) }    
    names(dist.lowest2) <- ok.chrs
    
    ## fix cases where two nested features were returned by choosing the lowest absolute distances for both features.
    res.nrst <- mapply(function(x,y) cbind(x,lowestDist=y),res.nrst,dist.lowest2, SIMPLIFY=F)
    res.nrst.i <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("res.nrst")) %dopar% {     
        mins <- tapply(abs(res.nrst[[x]][,"lowestDist"]),res.nrst[[x]][,"queryHits"],min); 
        res.nrst[[x]][mins[as.character(res.nrst[[x]][,"queryHits"])]==abs(res.nrst[[x]][,"lowestDist"]),]
    }
    names(res.nrst.i) <- ok.chrs
    
    if(!ismatrixFlag) {
        res.nrst.i <- sapply(ok.chrs,function(x) res.nrst.i[[x]][,"lowestDist"] )
    }
    
    return(res.nrst.i)
}

#' Generate a window size label.
#'
#' Function to generate asthetically pleasing window size label given an integer. This is one of the helper function used in \code{\link{getFeatureCounts}} & \code{\link{getFeatureCountsBig}}. 
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
	ind <- cut(abs(x), c(0, 1e3, 1e6, 1e9, 1e12), include.lowest = TRUE, right = FALSE, labels = FALSE)
	paste(x/c(1, 1e3, 1e6, 1e9, 1e12)[ind],c("bp", "Kb", "Mb", "Gb")[ind],sep="")
}

#' Resize a RangedData object.
#'
#' Function to resize a RangedData object by the given width, max space/chromosome size, and the boundary. This is one of the helper function used in \code{\link{getFeatureCountsBig}}. 
#'
#' @param rd a RangedData object
#' @param width the width of the resized ranges.
#' @param boundary same as fix parameter in \code{\link{resize}}. One of "start", "end", and "center". Default is "center".
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
#' alldata.rd <- makeRangedData(sites,soloStart=TRUE)
#' library(BSgenome.Hsapiens.UCSC.hg18)
#' resizeRangedData(alldata.rd,width=10000,spaceSizes=seqlengths(Hsapiens))
#' resizeRangedData(alldata.rd,width=10000,limitLess=T)
resizeRangedData <- function(rd,width=NULL,boundary="center",spaceSizes=NULL,spaceMin=1,limitLess=FALSE) {
	stopifnot(!is.null(width))
		
	new.rd <- resize(ranges(rd),width,fix=boundary)
	if(!as.logical(limitLess)) {
		stopifnot(!is.null(spaceSizes))
		new.rd <- restrict(new.rd, start=spaceMin)
		overEdge <- sapply(names(new.rd),function(x) which(end(new.rd[[x]])>spaceSizes[x]),simplify=FALSE)
		if (any(sapply(overEdge,length)>0)) {
			overEdge.chr <- names(which(sapply(overEdge,length)>0))
			for (f in overEdge.chr) {
				end(new.rd[[f]])[ as.numeric(unlist(overEdge[f]))] <- spaceSizes[f]
			} 
		}		
	}
	
	new.rd <- as(new.rd, "RangedData")
	values(new.rd) <- values(rd)
	return(new.rd)
}

#' Get counts of annotation within a defined window around each query range/position. 
#'
#' Given a query object and window size(s), the function finds all the rows in subject which are <= window size/2 distance away. If weights are assigned to each positions in the subject, then tallied counts are multiplied accordingly. If annotation object is large, spanning greater than 100 million rows, then \code{\link{getFeatureCountsBig}} is used which drops any weight column if specified or additional parameters passed to \code{\link{findOverlaps}}.
#'
#' @param sites.rd RangedData/GRanges object to be used as the query.
#' @param features.rd RangedData/GRanges object to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves as a prefix to windows sizes!
#' @param widths a named/numeric vector of window sizes to be used for casting a net around each position. Default: \code{c(1000,10000,1000000)}.
#' @param weightsColname if defined, weigh each row from features.rd when tallying up the counts.
#' @param doInChunks break up sites.rd into small pieces of chunkSize to perform the calculations. Default to FALSE. Useful if you are expecting to find great deal of overlap between sites.rd and features.rd.
#' @param chunkSize number of rows to use per chunk of sites.rd. Default to 10000. Only used if doInChunks=TRUE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to FALSE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param ... Additional parameters for \code{\link{findOverlaps}}. Default parameters passed are query, subject, select='all', maxgap=widths/2.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd. There will be a column for each width defined in widths parameter. If widths was a named vector i.e. c("100bp"=100,"1K"=1000), then the colname parameter will be pasted together with width name else default name will be generated by the function.
#'
#' @note 
#' \itemize{
#'   \item If the input sites.rd parameter is GRanges object, then it is converted to RangedData and then converted back to GRanges at the end since \code{\link{findOverlaps}} function operates much faster on RangedData objects. 
#'   \item Try not to use this function for >50 spaces unless you have tons fo memory. 
#'   \item If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#' }
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' head(sites)
#' alldata.rd <- makeRangedData(sites,soloStart=TRUE)
#' alldata.rd
#'
#' data(genes)
#' head(genes)
#' genes.rd <- makeRangedData(genes)
#' genes.rd
#'
#'
#' geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene")
#' geneCounts
#' # Parallel version of getFeatureCounts
#' geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene", parallel=T)
#' geneCounts
#' # For large annotations, use getFeatureCountsBig
getFeatureCounts <- function(sites.rd, features.rd, colnam=NULL,widths=c(1000,10000,1000000), weightsColname=NULL, doInChunks=FALSE, chunkSize=10000, parallel=FALSE, chromSizes=NULL, ...) {
	stopifnot(nrow(sites.rd)>0)
    stopifnot(nrow(features.rd)>0)
    
    grangesFlag <- FALSE
    if(class(sites.rd)=="GRanges" & class(features.rd)=="GRanges") {
    	grangesFlag <- TRUE
    	sites.rd <- as(sites.rd,"RangedData")
    	features.rd <- as(features.rd,"RangedData")
    }
    
    if (!any(names(sites.rd) %in% names(features.rd))) {
        stop("There are no spaces/chromosomes that are shared between the query (sites.rd) and subject (features.rd)")
    }
    if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
    }
    if(!is.null(chromSizes)) {
    	## for historical version...display a warning if this parameter is passed!
        warning("decrepit option: chromSizes parameter is no longer required and will be ignored!")
    }    
    
    if(nrow(features.rd)>=1e8) {
        message("Using getFeatureCountsBig() to get the job done!")
        getFeatureCountsBig(sites.rd, features.rd, colnam, widths)
    } else {
        if(doInChunks & chunkSize < nrow(sites.rd)) { # no need to execute all this if chunkSize is bigger than data size!!!           
            fun.call <- as.list(match.call())[-1]             
            extraParams <- fun.call[!names(fun.call) %in% c("sites.rd", "features.rd", "colnam", "widths", "weightsColname", "doInChunks", "chunkSize", "parallel")]
            total <- nrow(sites.rd)        
            starts <- seq(1,total,by=chunkSize)
            stops <- unique(c(seq(chunkSize ,total,by=chunkSize),total))
            stopifnot(length(starts)==length(stops))  
            message("Breaking up sites.rd into chunks of ",chunkSize)
            res <- RangedData()
            for(x in 1:length(starts)) {                
                if(length(extraParams)>1) {
                    res <- rbind(res,getFeatureCounts(sites.rd[starts[x]:stops[x],], features.rd, colnam, widths, weightsColname, parallel=parallel, eval(as.expression(extraParams))))
                } else {
                    res <- rbind(res,getFeatureCounts(sites.rd[starts[x]:stops[x],], features.rd, colnam, widths, weightsColname, parallel=parallel))
                }
            }
            
			if(grangesFlag) {
				re <- as(re,"GRanges")
			}
			
            return(res)
        } else {
            weighted <- ifelse(is.null(weightsColname),FALSE,TRUE)
            
            # only get labels if not supplied
            if(is.null(names(widths))) {
                names(widths) <- getWindowLabel(widths)
        	}
        	
            ## use only chromosomes that are present in both sites.rd and features.rd ##
            ok.chrs <- intersect(space(sites.rd),space(features.rd))
            features.rd <- features.rd[names(features.rd) %in% ok.chrs]
            
            query <- sites.rd[,-c(1:length(colnames(sites.rd)))]
            
            if(!parallel) { registerDoSEQ() }
            
            ## perform overlap analysis in parallel by windows
            allcounts <- foreach(x=iter(widths),.inorder=TRUE,.packages="IRanges",.export=c("query", "features.rd", "weighted", "weightsColname", "resizeRangedData")) %dopar% {            	
                
                if (weighted) {
                    res <- as.data.frame(as.matrix(findOverlaps(query, features.rd, select='all', maxgap=x/2, ...)))
                    res$weights <- features.rd[res$subjectHits,][[weightsColname]]
                    tapply(res$weights,res$queryHits,sum)            
                } else {
                    ## dont use countOverlaps() since it returns overlapping ranges from other spaces/chrs if it was a factor
                    as.table(findOverlaps(query, features.rd, select='all', maxgap=x/2, ...))  
                }
            }
            names(allcounts) <- names(widths)
            
            colnam <- cleanColname(colnam)
            
            ## add columns back to query object
            for (windowName in names(widths)) {
                columnName <- paste(colnam,names(widths[ windowName ] ),sep=".")
                sites.rd[[columnName]] <- 0
                sites.rd[[columnName]][as.numeric(names(allcounts[[windowName]]))] <- as.numeric(allcounts[[windowName]])
            }
            
            if(grangesFlag) {
				sites.rd <- as(sites.rd,"GRanges")
			}
			
            sites.rd
        }
    }
}

#' Clean the supplied string from punctuations and spaces.
#'
#' Function to clean the supplied string from punctuations and spaces so it can be used as column headings.
#'
#' @param x string or a vector to be cleaned.
#' @param description OPTIONAL string identifying the purpose supplied string in x to be displayed in the cleaning message.
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
cleanColname <- function(x, description="colnam") {
	if(any(grepl("[[:punct:]]",x) | grepl("[[:space:]]",x))) {
		message("Cleaning the supplied '",description,"'")
		x <- gsub("\\_+","_",gsub("[[:space:]]","_",gsub("[[:punct:]]","_",x)))
	}
	return(x)
}

#' Get counts of annotation within a defined window around each query range/position for large annotation objects spanning greater than 100 million rows. This is still in beta phase. 
#'
#' Given a query object and window size(s), the function finds all the rows in subject which are <= window size/2 distance away. This new range is then compared against the midpoints of subject to find any overlapping ranges and then tallied up. Note that here counting is done using midpoint of the ranges in subject instead of start-stop boundaries. 
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
getFeatureCountsBig <- function(sites.rd, features.rd, colnam=NULL, widths=c(1000,10000,1000000)) {
	stopifnot(nrow(sites.rd)>0)
    stopifnot(nrow(features.rd)>0)
    
    grangesFlag <- FALSE
    if(class(sites.rd)=="GRanges" & class(features.rd)=="GRanges") {
    	grangesFlag <- TRUE
    	sites.rd <- as(sites.rd,"RangedData")
    	features.rd <- as(features.rd,"RangedData")
    }
    
    if (!any(names(sites.rd) %in% names(features.rd))) {
        stop("There are no spaces/chromosomes that are shared between the query (sites.rd) and subject (features.rd)")
    }
    if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
    }

	# only get labels if not supplied
    if(is.null(names(widths))) {
        names(widths) <- getWindowLabel(widths)
    }
	
    ## use only chromosomes that are present in both sites.rd and features.rd ##
    ok.chrs <- intersect(space(sites.rd),space(features.rd))
    features.rd <- ranges(features.rd[names(features.rd) %in% ok.chrs])
    good.rows <- space(sites.rd) %in% ok.chrs
    
    colnam <- cleanColname(colnam)
    
    ## get counts of midpoints using findInterval and add columns back to query object
    for (windowName in names(widths)) {
        cat(".")
        columnName <- paste(colnam,names(widths[ windowName ] ),sep=".")
        sites.rd[[columnName]] <- 0
        
        query <- sites.rd[,-c(1:length(colnames(sites.rd)))]
        query <- ranges(resizeRangedData(query,width=widths[windowName],limitLess=TRUE))
        counts <- sapply(ok.chrs, function(x) {
                            features.mid <- sort((start(features.rd[[x]])+end(features.rd[[x]]))/2)                            
                            findInterval(end(query[[x]]),features.mid) - findInterval(start(query[[x]]),features.mid)
                        })             
        stopifnot(nrow(sites.rd)==sum(unlist(lapply(counts,length)))) 
        sites.rd[[columnName]][good.rows] <- unsplit(counts, space(sites.rd)[good.rows], drop = TRUE)
    }
    
    if(grangesFlag) {
		sites.rd <- as(sites.rd,"GRanges")
	}
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
#' @param strand.colnam column name from features.rd to be used for retrieving the feature orientation. By default this is NULL assuming that features.rd has a column that includes the word 'strand' somewhere in it. If it doesn't the function will assume the supplied annotation has no orientation information '*'. Not required if asBool=TRUE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to FALSE. Not applicable when asBool=T. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param ... Additional parameters for \code{\link{findOverlaps}}.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd.
#'
#' @note 
#' \itemize{
#'   \item If the input sites.rd parameter is GRanges object, then it is converted to RangedData and then converted back to GRanges at the end since \code{\link{findOverlaps}} function operates much faster on RangedData objects. 
#'   \item Try not to use this function for >50 spaces unless you have tons fo memory. 
#'   \item If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#' }
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}}, \code{\link{getNearestFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' head(sites)
#' alldata.rd <- makeRangedData(sites,soloStart=TRUE)
#' alldata.rd
#'
#' data(genes)
#' head(genes)
#' genes.rd <- makeRangedData(genes)
#' genes.rd
#'
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene")
#' InGenes
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene",asBool=TRUE)
#' InGenes
#' # Parallel version of getSitesInFeature
#' InGenes <- getSitesInFeature(alldata.rd,genes.rd,"InGene",asBool=TRUE,parallel=TRUE)
#' InGenes
getSitesInFeature <- function(sites.rd, features.rd, colnam=NULL, asBool=F, feature.colnam=NULL, strand.colnam=NULL, parallel=FALSE, ...) {    
	stopifnot(nrow(sites.rd)>0)
    stopifnot(nrow(features.rd)>0)
    
    grangesFlag <- FALSE
    if(class(sites.rd)=="GRanges" & class(features.rd)=="GRanges") {
    	grangesFlag <- TRUE
    	sites.rd <- as(sites.rd,"RangedData")
    	features.rd <- as(features.rd,"RangedData")
    }
    
    if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
    }
    
    if (!any(names(sites.rd) %in% names(features.rd))) {
        stop("There are no spaces/chromosomes that are shared between the query (sites.rd) and subject (features.rd)")
    }
    
    if (asBool) {
        res <- sites.rd %in% features.rd
        sites.rd[[colnam]] <- unlist(res)
    } else {
        ## find overlaps!
        res <- as.data.frame(as.matrix(findOverlaps(sites.rd,features.rd,...)))
        stopifnot(!any(is.na(res)))                    
        
        ## get feature names and strand for identifying overlapping genes
        if(is.null(feature.colnam)) {
            featureName <- getRelevantCol(colnames(features.rd),c("name","featureName"),"featureName",multiple.ok=TRUE)
            feature.colnam <- colnames(features.rd)[featureName][1]
        }
        
        skipStrandCalc <- FALSE
        if(is.null(strand.colnam)) {
            strandCol <- grep('strand',colnames(features.rd),ignore.case=TRUE,value=TRUE)[1]
            if(is.na(strandCol)) {
                skipStrandCalc <- TRUE
            } else {
                strand.colnam <- strandCol
            }
        }
                
        res$features <- unlist(values(features.rd))[feature.colnam][,1][res$subjectHits]
        if(skipStrandCalc) {
            res$strands <- "*"
        } else {
            res$strands <- unlist(values(features.rd))[strand.colnam][,1][res$subjectHits]
        }
        
        ## collapse rows where query returned two hits with the same featureNames due to alternative splicing or something else.
        res <- unique(res[,c("queryHits","features","strands")]) 
        
        counts <- table(res$query)    
        res$multipleHit <- res$queryHits %in% as.numeric(names(counts[counts>1]))
        res$featureName <- res$features
        res$strand <- res$strands
        
        if(!parallel) { registerDoSEQ() }
        
        if(any(res$multipleHit)) {
            res.unique <- unique(subset(res,!multipleHit,select=c("queryHits","featureName","strand"),drop=TRUE))
            res <- subset(res,multipleHit,drop=TRUE)        
            
            combos <- split(res$features,res$queryHits)
            tmp <- foreach(x=iter(combos),.inorder=TRUE) %dopar% { paste(unique(x),collapse=",") }
            names(tmp) <- names(combos)        
            res[,"featureName"] <- unlist(tmp[as.character(res$queryHits)])
            
            if(skipStrandCalc) {
                res[,"strand"] <- "*"
            } else {
                combos <- split(res$strands,res$queryHits)
                tmp <- foreach(x=iter(combos),.inorder=TRUE) %dopar% { paste(unique(x),collapse=",") }
                names(tmp) <- names(combos)        
                res[,"strand"] <- unlist(tmp[as.character(res$queryHits)])
            }
            
            res <- unique(res[,c("queryHits","featureName","strand")])
            res <- rbind(res,res.unique)
        } else {
            res <- unique(subset(res,!multipleHit,select=c("queryHits","featureName","strand"),drop=TRUE))
        }
        
        ## for collapsed rows take the uniques and add them back to sites.rd        
        colnam <- cleanColname(colnam)
        
        sites.rd[[colnam]] <- FALSE
        sites.rd[[paste(colnam,"Ort",sep="")]] <- NA
        sites.rd[[colnam]][res$query] <- as.character(res$featureName)
        sites.rd[[paste(colnam,"Ort",sep="")]][res$query] <- as.character(res$strand)
    }
    
    if(grangesFlag) {
		sites.rd <- as(sites.rd,"GRanges")
	}
	
    sites.rd
}

#' Annotate a RangedData/GRanges object using one of annotation functions. 
#'
#' This is a wrapper function which calls one of following functions depending on annotType parameter: \code{\link{getFeatureCounts}}, \code{\link{getNearestFeature}}, code{\link{getSitesInFeature}} 
#'
#' @param annotType one of following: within, nearest, counts.
#' @param ... Additional parameters to be passed to the respective annotation function.
#'
#' @return a RangedData/GRanges object with new annotation columns appended at the end of sites.rd.
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{makeGRanges}}, \code{\link{getFeatureCounts}}, \code{\link{getNearestFeature}}, code{\link{getSitesInFeature}}.
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to RangedData/GRanges object
#' data(sites)
#' head(sites)
#' alldata.rd <- makeRangedData(sites,soloStart=TRUE)
#' alldata.rd
#'
#' data(genes)
#' head(genes)
#' genes.rd <- makeRangedData(genes)
#' genes.rd
#'
#' library(BSgenome.Hsapiens.UCSC.hg18)
#'
#' doAnnotation(annotType="within",alldata.rd,genes.rd,"InGene")
#' doAnnotation(annotType="counts",alldata.rd,genes.rd,"NumOfGene")
#' doAnnotation(annotType="nearest",alldata.rd,genes.rd,"NearestGene")
doAnnotation <- function(annotType=NULL,...) {    
    if(is.null(annotType)) {
        stop("Please define the annotType parameter to identify which type of annotation to perform: within, nearest, counts")
    }
	
	switch(EXPR = annotType,
		within = getSitesInFeature(...),
		nearest = getNearestFeature(...),
		counts = getFeatureCounts(...),
		stop("Invalid annoType")
	)
}
