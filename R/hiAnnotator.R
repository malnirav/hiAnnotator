#' Initiate UCSC genome browser session given the freeze argument.
#'
#' @param freeze one of following: hg18, mm8, rheM, etc. Default is hg18.
#'
#' @return browser session object compatible with rtracklayer functions.
#'
#' @seealso \code{\link{getUCSCtable}}, \code{\link{makeRangedData}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @examples
#' session<-makeUCSCsession()
#' genome(session)
#' session<-makeUCSCsession("mm8")
#' genome(session)
makeUCSCsession<-function(freeze="hg18") {
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
#' @examples
#' refflat<-getUCSCtable("refFlat","RefSeq Genes") ## same as session<-makeUCSCsession(); refflat<-getUCSCtable("refFlat", "RefSeq Genes", bsession=session, freeze="hg18")
#' head(refflat)
getUCSCtable<-function(tableName,trackName,bsession=NULL,freeze="hg18",...) {
    if(is.null(bsession)) { 
        bsession<-makeUCSCsession(freeze)
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
#' @seealso \code{\link{makeRangedData}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @examples
#' data(sites)
#' names(sites)
#' getRelevantCol(names(sites),c("chr","chromosome","tname","space","chrom","contig"),"space")
#' getRelevantCol(names(sites),c("ort","orientation","strand"),"strand")
getRelevantCol<-function(col.names,col.options,col.type,multiple.ok=FALSE) {
    answer<-unique(as.numeric(unlist(sapply(col.options,function(x) grep(x,col.names,ignore.case=TRUE)))))
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
#' @param x data frame to be converted into RangedData object
#' @param positionsOnly boolean flag indicating to return only position based data or everything from the data frame. Defaults to FALSE.
#' @param soloStart flag denoting whether only one position based column is available. In other words, only starts are present and no stops. Default=FALSE.
#' @param chromCol use the defined column name for space/chromosome based data from the data frame. Defaults to NULL.
#' @param strandCol use the defined column name for strand or orientation from the data frame. Defaults to NULL.
#' @param startCol use the defined column name for start coordinate from the data frame. Defaults to NULL.
#' @param stopCol use the defined column name for stop coordinate from the data frame. Defaults to NULL and not required if soloStart=TRUE.
#'
#' @return a RangedData object converted from x.
#'
#' @seealso \code{\link{getNearestFeature}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}.
#'
#' @examples
#' # Convert a data frame to RangedData object
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
makeRangedData<-function(x,positionsOnly=FALSE,soloStart=FALSE,chromCol=NULL,strandCol=NULL,startCol=NULL,stopCol=NULL) {
    ## set column names for space and strand if not provided ##
    if(is.null(chromCol)) {
        colIndex<-getRelevantCol(names(x),c("chr","chromosome","tname","space","chrom","contig"),"space")
        names(x)[colIndex]<-"space"
    } else {
        names(x)[names(x)==chromCol]<-"space"
    }
    
    if(is.null(strandCol)) {
        colIndex<-getRelevantCol(names(x),c("ort","orientation","strand"),"strand")
        names(x)[colIndex]<-"strand"
    } else {
        names(x)[names(x)==strandCol]<-"strand"
    }    
        
    if(is.null(startCol)) {
        startCol<-getRelevantCol(names(x),c("position","intsite","txstart","start","chromstart"),"start",multiple.ok=TRUE)
        names(x)[startCol[1]]<-"start" 
    } else {
        names(x)[names(x)==startCol]<-"start"
    }
    
    ## only do stop if soloStart=F ##
    if(!as.logical(soloStart) & is.null(stopCol)) {
        stopCol<-getRelevantCol(names(x),c("txend","end","stop","chromend"),"end",multiple.ok=TRUE)
        names(x)[stopCol[1]]<-"end" 
    } else {
        names(x)[names(x)==stopCol]<-"end"
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
        for (y in names(which(factorCols))) { x[,y]<-as.character(x[,y]); if(!any(is.na(suppressWarnings(as.numeric(x[,y]))))) { x[,y]<-as.numeric(x[,y]) } }        
    }

    ## if start and end coordinates are present, sort by midpoint ##
    ## else if only single coordinate is present, then add the end column and sort ##
    if(length(startCol)>0 & length(stopCol)>0) {       
        if(any(is.na(x$end))) {
            stop("NAs found in column containing end positions")
        }
        x$mid<-with(x,(start+end)/2)
        x<-orderBy(~space+mid,x)
        x$mid<-NULL
    } else {  
        x<-orderBy(~space+start,x)
        x$end<-x$start
    }
    
    if(as.logical(positionsOnly)) {
        sites.rd<-as(x[,c("space","start","end","strand")],"RangedData")
    } else {
        sites.rd<-as(x,"RangedData")
    }
    
    sites.rd
}

#' Get nearest annotation boundary for a position range. 
#'
#' Given a query object, the function retrieves the nearest feature and its properties from a subject and then appends them as new columns within the query object. When used in genomic context, the function can be used to retrieve a nearest gene 5' or 3' end relative to genomic position of interest.
#'
#' @param sites.rd RangedData object to be used as the query.
#' @param features.rd RangedData obj to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves a core!
#' @param side boundary of annotation to use to calculate the nearest distance. Options are '5p','3p', or the default 'either'.
#' @param feature.colnam column name from features.rd to be used for retrieving the nearest feature name. By default this is NULL assuming that features.rd has a column that includes the word 'name' somewhere in it.
#' @param strand.colnam column name from features.rd to be used for retrieving the nearest feature's orientation. By default this is NULL assuming that features.rd has a column that includes the word 'strand' somewhere in it. If it doesn't the function will assume the supplied annotation is in '+' orientation (5' -> 3'). The same applies to strand column in sites.rd.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to FALSE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#'
#' @return a RangedData object with new annotation columns appended at the end of sites.rd.
#'
#' @note Try not to use this function for >50 spaces unless you have tons fo memory. If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}, \code{\link{get2NearestFeature}}.
#'
#' @examples
#' # Convert a data frame to RangedData object
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
    if (!any(names(sites.rd) %in% names(features.rd))) {
        stop("There are no spaces/chromosomes that are shared between the query (sites.rd) and subject (features.rd)")
    }
    if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
    }
    ## use only chromosomes that are present in both sites.rd and features.rd ##
    ok.chrs <- intersect(space(sites.rd),space(features.rd))
    features.rd$usablerows<-space(features.rd) %in% ok.chrs
    features.rd <- subset(features.rd,usablerows,drop=TRUE)
    features.rd$usablerows<-NULL
    
    if(is.null(feature.colnam)) {
        featureName<-getRelevantCol(colnames(features.rd),c("name","featureName"),"featureName",multiple.ok=TRUE)
        feature.colnam<-colnames(features.rd)[featureName][1]  
    }
    
    ## do a check of strand column in sites.rd ##
    answer <- unique(as.numeric(unlist(sapply(c("strand","orientation","ort"), function(x) grep(x,colnames(sites.rd), ignore.case = TRUE)))))
    strandCol<-colnames(sites.rd)[answer][1]
    if(is.na(strandCol)) { 
        message("No orientation column found in sites.rd. Using '+' as default.")
        sites.rd$strand<-"+" 
    }
    
    ## do a check of strand column in features.rd ##
    answer <- unique(as.numeric(unlist(sapply(c("strand","orientation","ort"), function(x) grep(x,colnames(features.rd), ignore.case = TRUE)))))
    strandCol<-colnames(features.rd)[answer][1]
    if(is.null(strand.colnam) & is.na(strandCol)) { 
        message("No orientation column found in features.rd. Using '+' as default.")
        features.rd$strand<-"+"
    }
    if(!is.na(strandCol) & strandCol!="strand") { features.rd$strand<-unlist(values(features.rd))[strandCol] }
    
    ## convert any factor columns to character to avoid downstream issues with NAs and unlisting of CompressedCharacterList object ##
    factorCols <- sapply(colnames(features.rd),function(x) class(features.rd[[x]]))=="factor"
    if(any(factorCols)) {
        for (x in names(which(factorCols))) { features.rd[[x]]<-as.character(features.rd[[x]]) }        
    }
    
    vals.s <- values(features.rd)
    
    query <- ranges(sites.rd)

    if (side=='5p')
        subject <- flank(ranges(features.rd),start=split(features.rd$strand=="+",features.rd$space),width=-1) ##get only 5 prime sides of features
    
    if (side=='3p')
        subject <- flank(ranges(features.rd),start=split(features.rd$strand=="-",features.rd$space),width=-1) ##get only 3 prime sides of features
    
    if (side=='either')
        subject <- ranges(features.rd)
    
    if(!parallel) { registerDoSEQ() }
    
    # first get the nearest indices
    res <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject"),.packages="IRanges") %dopar% { as.matrix(nearest(query[[x]], subject[[x]], select="all")) }     
    names(res)<-ok.chrs
    stopifnot(identical(lapply(query,length)[ok.chrs],lapply(res,function(x) length(unique(x[,"query"])) )[ok.chrs])) ## check for safety
    
    # check if >1 nearest matches found, get the indices of the feature with shortest distance to 5p/3p
    res <- getLowestDists(query,subject,subjectOrt=vals.s[,"strand"],ok.chrs=ok.chrs,res.nrst=res,side=side,cores.use=1)
    
    # fix any cases where matrix got converted to integer due to only value found
    if(any(unlist(lapply(res,class))!="matrix")) {
        tofix <- names(which(unlist(lapply(res,class))!="matrix"))
        for (i in tofix) { res[[i]] <- t(res[[i]]) }
    }
    
    # for the feature of shortest indices, get the names, and strand attributes
    featureName <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("vals.s","res","feature.colnam")) %dopar% { vals.s[[x]][res[[x]][,"subject"],feature.colnam] }     
    
    ort <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("vals.s","res","feature.colnam")) %dopar% { vals.s[[x]][res[[x]][,"subject"],"strand"] }     
    
    names(featureName) <- names(ort) <- ok.chrs
    
    stopifnot(identical(lapply(res,nrow),lapply(featureName,length))) ## check for safety
    stopifnot(identical(lapply(res,nrow),lapply(ort,length))) ## check for safety
    
    res.i <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("ort","res","featureName")) %dopar% {
        toprune <- as.data.frame(res[[x]])
        toprune$featureName <- featureName[[x]]
        toprune$ort <- ort[[x]]
        toprune <- unique(toprune[,-2])
        counts <- table(toprune$query)
        ismulti <- toprune$query %in% as.numeric(names(counts[counts>1]))
        if(any(ismulti)) {
            goods <- toprune[!ismulti,]
            toprune <- toprune[ismulti,]
            tempStore <- with(toprune,tapply(featureName,query,paste,collapse=","))
            toprune$featureName <- as.character(tempStore[as.character(toprune$query)])            
            tempStore <- with(toprune,tapply(ort,query,paste,collapse=","))
            toprune$ort <- as.character(tempStore[as.character(toprune$query)])            
            tempStore <- with(toprune,sapply(tapply(lowestDist,query,abs),min)) ## if a site falls exactly between two genes pick one the abs(lowest)
            toprune$lowestDist <- as.numeric(tempStore[as.character(toprune$query)])
            toprune <- rbind(goods,unique(toprune))
        }
        return(orderBy(~query,toprune))
    }  
    names(res.i) <- ok.chrs
    
    prefix <- ifelse(side=="either","",side)
    good.rows <- space(sites.rd) %in% ok.chrs
    
    sites.rd[[paste(prefix,colnam,sep="")]][good.rows] <- unsplit(lapply(res.i,"[[","featureName"), space(sites.rd)[good.rows], drop = TRUE)
    sites.rd[[paste(prefix,colnam,"Ort",sep="")]][good.rows] <- unsplit(lapply(res.i,"[[","ort"), space(sites.rd)[good.rows], drop = TRUE)
    sites.rd[[paste(prefix,colnam,"Dist",sep="")]][good.rows] <- unsplit(lapply(res.i,"[[","lowestDist"), space(sites.rd)[good.rows], drop = TRUE)    
    sites.rd
}

#' Get two nearest upstream and downstream annotation boundary for a position range. 
#'
#' Given a query object, the function retrieves the two nearest feature upstream and downstream along with their properties from a subject and then appends them as new columns within the query object. When used in genomic context, the function can be used to retrieve two nearest gene upstream and downstream of the genomic position of interest.
#'
#' @param sites.rd RangedData object to be used as the query.
#' @param features.rd RangedData obj to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves a core!
#' @param side boundary of annotation to use to calculate the nearest distance. Options are '5p','3p', or the default 'either'.
#' @param feature.colnam column name from features.rd to be used for retrieving the nearest feature name. By default this is NULL assuming that features.rd has a column that includes the word 'name' somewhere in it.
#' @param strand.colnam column name from features.rd to be used for retrieving the nearest feature's orientation. By default this is NULL assuming that features.rd has a column that includes the word 'strand' somewhere in it. If it doesn't the function will assume the supplied annotation is in '+' orientation (5' -> 3'). The same applies to strand column in sites.rd.
#'
#' @return a RangedData object with new annotation columns appended at the end of sites.rd.
#'
#' @note For cases where a position is at the edge and there are no feature up/down stream since it would fall off the chromosome, the function simply returns the nearest feature. In addition, if there are multiple locations where a query falls into, the function arbitrarily chooses one to serve as the nearest feature, then reports 2 upstream & downstream feature. That may occasionally yield features which are the same upstream and dowstream, which is commonly encountered when studying spliced genes or phenomena related to it. 
#'
#' @seealso \code{\link{getNearestFeature}}, \code{\link{makeRangedData}}, \code{\link{getFeatureCounts}}, \code{\link{getSitesInFeature}}.
#'
#' @examples
#' # Convert a data frame to RangedData object
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
    if (!any(names(sites.rd) %in% names(features.rd))) {
        stop("There are no spaces/chromosomes that are shared between the query (sites.rd) and subject (features.rd)")
    }
    if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
    }
    ## use only chromosomes that are present in both sites.rd and features.rd ##
    ok.chrs <- intersect(space(sites.rd),space(features.rd))
    features.rd$usablerows<-space(features.rd) %in% ok.chrs
    features.rd <- subset(features.rd,usablerows,drop=TRUE)
    features.rd$usablerows<-NULL
    
    if(is.null(feature.colnam)) {
        featureName<-getRelevantCol(colnames(features.rd),c("name","featureName"),"featureName",multiple.ok=TRUE)
        feature.colnam<-colnames(features.rd)[featureName][1]  
    }
    
    ## do a check of strand column in sites.rd ##
    answer <- unique(as.numeric(unlist(sapply(c("strand","orientation","ort"), function(x) grep(x,colnames(sites.rd), ignore.case = TRUE)))))
    strandCol<-colnames(sites.rd)[answer][1]
    if(is.na(strandCol)) { 
        message("No orientation column found in sites.rd. Using '+' as default.")
        sites.rd$strand<-"+" 
    }
    
    ## do a check of strand column in features.rd ##
    answer <- unique(as.numeric(unlist(sapply(c("strand","orientation","ort"), function(x) grep(x,colnames(features.rd), ignore.case = TRUE)))))
    strandCol<-colnames(features.rd)[answer][1]
    if(is.null(strand.colnam) & is.na(strandCol)) { 
        message("No orientation column found in features.rd. Using '+' as default.")
        features.rd$strand<-"+"
    }
    if(!is.na(strandCol) & strandCol!="strand") { features.rd$strand<-unlist(values(features.rd))[strandCol] }

    ## convert any factor columns to character to avoid downstream issues with NAs and unlisting of CompressedCharacterList object ##
    factorCols <- sapply(colnames(features.rd),function(x) class(features.rd[[x]]))=="factor"
    if(any(factorCols)) {
        for (x in names(which(factorCols))) { features.rd[[x]]<-as.character(features.rd[[x]]) }        
    }
    
    vals.q <- values(sites.rd) 
    vals.s <- values(features.rd)
    
    query <- ranges(sites.rd)

    if (side=='5p')
        subject <- flank(ranges(features.rd),start=split(features.rd$strand=="+",features.rd$space),width=-1) ##get only 5 prime sides of features
    
    if (side=='3p')
        subject <- flank(ranges(features.rd),start=split(features.rd$strand=="-",features.rd$space),width=-1) ##get only 3 prime sides of features
    
    if (side=='either')
        subject <- ranges(features.rd)
    
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
    
    u1 <- sapply(ok.chrs, function(x) ifelse(dist.nrst[[x]]<0, res[[x]], ifelse(vals.q[[x]][,"strand"]=="+",res.left1[[x]],res.right1[[x]])))   
    d1 <- sapply(ok.chrs, function(x) ifelse(dist.nrst[[x]]<0, ifelse(vals.q[[x]][,"strand"]=="+",res.right1[[x]],res.left1[[x]]), res[[x]]))
    u2 <- sapply(ok.chrs, function(x) ifelse(dist.nrst[[x]]<0, ifelse(vals.q[[x]][,"strand"]=="+",res.left1[[x]],res.right1[[x]]), ifelse(vals.q[[x]][,"strand"]=="+",res.left2[[x]],res.right2[[x]])))
    d2 <- sapply(ok.chrs, function(x) ifelse(dist.nrst[[x]]<0, ifelse(vals.q[[x]][,"strand"]=="+",res.right2[[x]],res.left2[[x]]), ifelse(vals.q[[x]][,"strand"]=="+",res.right1[[x]],res.left1[[x]])))

    prefix <- ifelse(side=="either","Either",side)
    good.rows <- space(sites.rd) %in% ok.chrs
    
    message("u = upstream, d = downstream")
    message("thinking concept: u2.....u1.....intSite(+).....d1.....d2")
    message("thinking concept: d2.....d1.....intSite(-).....u1.....u2")
    
    for (f in c("u1","u2","d1","d2")) {
        message(f)
        res.nrst<-get(f)
        res.nrst <- sapply(ok.chrs, function(x) ifelse(res.nrst[[x]]>length(subject[[x]]) | res.nrst[[x]]<1,res[[x]],res.nrst[[x]])) 
        features <- sapply(ok.chrs, function(x) vals.s[[x]][res.nrst[[x]],feature.colnam])
        ort <- sapply(ok.chrs, function(x) vals.s[[x]][res.nrst[[x]],"strand"]) 
        dists <- getLowestDists(query,subject,vals.s[,"strand"],ok.chrs,res.nrst,side,2)
        
        if (f == "u1") { coldef<-paste(prefix,colnam,"upStream1",sep=".") }
        if (f == "u2") { coldef<-paste(prefix,colnam,"upStream2",sep=".") }
        if (f == "d1") { coldef<-paste(prefix,colnam,"downStream1",sep=".") }
        if (f == "d2") { coldef<-paste(prefix,colnam,"downStream2",sep=".") }
        
        sites.rd[[coldef]][good.rows] <- unsplit(features, space(sites.rd)[good.rows], drop = TRUE)
        sites.rd[[paste(coldef,"Ort",sep=".")]][good.rows] <- unsplit(ort, space(sites.rd)[good.rows], drop = TRUE)
        sites.rd[[paste(coldef,"Dist",sep=".")]][good.rows] <- unsplit(dists, space(sites.rd)[good.rows], drop = TRUE)    
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
getLowestDists<-function(query=NULL,subject=NULL,subjectOrt=NULL,ok.chrs=NULL,res.nrst=NULL,side="either",cores.use=1) {
    if(is.null(query) | is.null(subjectOrt) | is.null(subject) | is.null(ok.chrs) | is.null(res.nrst)) {
        stop("One of following is null: query, subjectOrt, subject, ok.chrs, res.nrst")
    }    
    
    if(any(!names(subject) %in% names(subjectOrt))) {
        stop("There are missing spaces/chromosomes in subjectOrt which are present in subject")
    }
    
    ## convert res.nrst to a matrix if its a list of integers/numeric resulting from using the default nearest() used in get2NearestFeature()
    ismatrixFlag <- TRUE
    if(all(unlist(lapply(res.nrst,class))!="matrix")) {
        res.nrst <- lapply(res.nrst,function(x) as.matrix(cbind(query=1:length(x),subject=x)) )
        ismatrixFlag <- FALSE
    }

    if(any(unlist(sapply(ok.chrs,function(x) res.nrst[[x]][,"subject"]>length(subject[[x]]) | res.nrst[[x]][,"subject"]<1)))) {
        stop("Subject indicies out of vector range in res.nrst ... this can introduce NAs and cause science to fail!")
    }
    
    if(cores.use==1) { registerDoSEQ() }
    
    if (side=="either") {
        ## get the lowest dist to either annot boundary from 5p side of the query
        dist.s <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { start(query[[x]][res.nrst[[x]][,"query"]])-start(subject[[x]][res.nrst[[x]][,"subject"]]) }     
        dist.e <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { start(query[[x]][res.nrst[[x]][,"query"]])-end(subject[[x]][res.nrst[[x]][,"subject"]]) }     
        names(dist.s) <- names(dist.e) <- ok.chrs

        dist5p <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("dist.s","dist.e")) %dopar% { ifelse(abs(dist.s[[x]])<abs(dist.e[[x]]),dist.s[[x]],dist.e[[x]]) }

        ## get the lowest dist to either annot boundary from 3p side of the query
        dist.s <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { end(query[[x]][res.nrst[[x]][,"query"]])-start(subject[[x]][res.nrst[[x]][,"subject"]]) }     
        dist.e <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { end(query[[x]][res.nrst[[x]][,"query"]])-end(subject[[x]][res.nrst[[x]][,"subject"]]) }     
        names(dist.s) <- names(dist.e) <- ok.chrs

        dist3p <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("dist.s","dist.e")) %dopar% { ifelse(abs(dist.s[[x]])<abs(dist.e[[x]]),dist.s[[x]],dist.e[[x]]) }               
    } else {
        ## get the lowest dist to annot boundary from 5p side of the query
        dist5p <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { start(query[[x]][res.nrst[[x]][,"query"]])-start(subject[[x]][res.nrst[[x]][,"subject"]]) }
        
        ## get the lowest dist to annot boundary from 3p side of the query
        dist3p <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("query","subject","res.nrst"),.packages="IRanges") %dopar% { end(query[[x]][res.nrst[[x]][,"query"]])-start(subject[[x]][res.nrst[[x]][,"subject"]]) }
    }
    
    names(dist3p) <- names(dist5p) <- ok.chrs
        
    ## get the lowest distance from the lowest 5p or 3p of the query!
    dist.lowest <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("dist5p","dist3p")) %dopar% { ifelse(abs(dist5p[[x]])<abs(dist3p[[x]]),dist5p[[x]],dist3p[[x]]) }
    names(dist.lowest) <- ok.chrs
    
    stopifnot(identical(unlist(lapply(res.nrst[ok.chrs],nrow)),unlist(lapply(dist.lowest,length)))) ## check for safety
    
    ## fix signs to match biological upstream or downstream
    dist.lowest2 <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("subjectOrt","res.nrst","dist.lowest")) %dopar% { ifelse(subjectOrt[[x]][res.nrst[[x]][,"subject"]]=="-",-dist.lowest[[x]],dist.lowest[[x]]) }    
    names(dist.lowest2) <- ok.chrs
    
    res.nrst <- mapply(function(x,y) cbind(x,lowestDist=y),res.nrst,dist.lowest2, SIMPLIFY=F)
    res.nrst.i <- foreach(x=iter(ok.chrs),.inorder=TRUE,.export=c("res.nrst")) %dopar% {     
        mins <- tapply(abs(res.nrst[[x]][,"lowestDist"]),res.nrst[[x]][,"query"],min); 
        res.nrst[[x]][mins[as.character(res.nrst[[x]][,"query"])]==abs(res.nrst[[x]][,"lowestDist"]),]
    }
    names(res.nrst.i)<-ok.chrs
    
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
#' @seealso \code{\link{getFeatureCounts}}, \code{\link{makeRangedData}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @examples
#' # Convert a data frame to RangedData object
#' getWindowLabel(c(0,1e7,1e3,1e6,2e9))
getWindowLabel<-function(x) {
	ind <- cut(abs(x), c(0, 1e3, 1e6, 1e9, 1e12), include.lowest = TRUE, right = FALSE, labels = FALSE)
	paste(x/c(1, 1e3, 1e6, 1e9, 1e12)[ind],c("bp", "Kb", "Mb", "Gb")[ind],sep="")
}

#' Resize a RangedData object.
#'
#' Function to resize a RangedData object by the given width, max space/chromosome size, and the boundary. This is one of the helper function used in \code{\link{getFeatureCounts}} & \code{\link{getFeatureCountsBig}}. 
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
#' @seealso \code{\link{getFeatureCounts}}, \code{\link{makeRangedData}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @examples
#' # Convert a data frame to RangedData object
#' library(BSgenome.Hsapiens.UCSC.hg18)
#' data(sites)
#' alldata.rd <- makeRangedData(sites,soloStart=TRUE)
#' resizeRangedData(alldata.rd,width=10000,spaceSizes=seqlengths(Hsapiens))
#' resizeRangedData(alldata.rd,width=10000,limitLess=T)
resizeRangedData<-function(rd,width=NULL,boundary="center",spaceSizes=NULL,spaceMin=1,limitLess=FALSE) {
	stopifnot(!is.null(width))
		
	new.rd <- resize(ranges(rd),width,fix=boundary)
	if(!as.logical(limitLess)) {
		stopifnot(!is.null(spaceSizes))
		new.rd <- restrict(new.rd, start=spaceMin)
		overEdge<-sapply(names(new.rd),function(x) which(end(new.rd[[x]])>spaceSizes[x]),simplify=FALSE)
		if (any(sapply(overEdge,length)>0)) {
			overEdge.chr<-names(which(sapply(overEdge,length)>0))
			for (f in overEdge.chr) {
				end(new.rd[[f]])[ as.numeric(unlist(overEdge[f]))]<-spaceSizes[f]
			} 
		}		
	}
	
	new.rd <- as(new.rd, "RangedData")
	values(new.rd) <- values(rd)
	return(new.rd)
}

#' Get counts of annotation within a defined window around each query range/position. 
#'
#' Given a query object and window size(s), the function first augments the ranges within the query by flanking starts and stops with window width. Therefore, a start of 12 and end of 14 with width 10 will yield a range of 8,17. This new range is then compared against the subject to find any overlapping ranges and then tallied up. If weights are assigned to each positions in the subject, then tallied counts are multiplied accordingly. If annotation object is large, spanning greater than 100 million rows, then getFeatureCountsBig is used which drops any weight column if specified or additional parameters passed to \code{\link{findOverlaps}}.
#'
#' @param sites.rd RangedData object to be used as the query.
#' @param features.rd RangedData obj to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves as a prefix to windows sizes!
#' @param chromSizes named vector of chromosome/space sizes to be used for testing if a position is off the mappable region.
#' @param widths a named/numeric vector of window sizes to be used for casting a net around each position. Default: \code{c(1000,10000,1000000)}
#' @param weightsColname if defined, weigh each row from features.rd when performing the counts
#' @param doInChunks break up sites.rd into small pieces of chunkSize to perform the calculations. Default to FALSE. Useful if you are expecting to find great deal of overlap between sites.rd and features.rd.
#' @param chunkSize number of rows to use per chunk of sites.rd. Default to 10000. Only used if doInChunks=TRUE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to FALSE. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param ... Additional parameters for \code{\link{findOverlaps}}.
#'
#' @return a RangedData object with new annotation columns appended at the end of sites.rd. There will be a column for each width defined in widths parameter. If widths was a named vector i.e. c("100bp"=100,"1K"=1000), then the colname parameter will be pasted together with width name else default name will be generated by the function.
#'
#' @note Try not to use this function for >50 spaces unless you have tons fo memory. If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}.
#'
#' @examples
#' # Convert a data frame to RangedData object
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
#' geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene",seqlengths(Hsapiens))
#' geneCounts
#' # Parallel version of getFeatureCounts
#' geneCounts <- getFeatureCounts(alldata.rd,genes.rd,"NumOfGene",seqlengths(Hsapiens), parallel=T)
#' geneCounts
#' # For large annotations, use getFeatureCountsBig
getFeatureCounts <- function(sites.rd, features.rd, colnam=NULL, chromSizes=NULL, widths=c(1000,10000,1000000), weightsColname=NULL, doInChunks=FALSE, chunkSize=10000, parallel=FALSE, ...) {
    if (!any(names(sites.rd) %in% names(features.rd))) {
        stop("There are no spaces/chromosomes that are shared between the query (sites.rd) and subject (features.rd)")
    }
    if(is.null(colnam)) {
        stop("Please define the colnam parameter for the new column(s) to be appended.")
    }
    if(is.null(chromSizes)) {
        stop("Please provide chromosome sizes. chromSizes=c('chr1'=247249719, 'chr2'=242951149...)")
    }
    
    if(nrow(features.rd)>=1e8) {
        message("Using getFeatureCountsBig() to get the job done!")
        getFeatureCountsBig(sites.rd, features.rd, colnam, widths)
    } else {
        if(doInChunks & chunkSize < nrow(sites.rd)) { # no need to execute all this if chunkSize is bigger than data size!!!           
            fun.call <- as.list(match.call())[-1]             
            extraParams <- fun.call[!names(fun.call) %in% c("sites.rd", "features.rd", "colnam", "chromSizes", "widths", "weightsColname", "doInChunks", "chunkSize", "parallel")]
            total <- nrow(sites.rd)        
            starts <- seq(1,total,by=chunkSize)
            stops <- unique(c(seq(chunkSize ,total,by=chunkSize),total))
            stopifnot(length(starts)==length(stops))  
            message("Breaking up sites.rd into chunks of ",chunkSize)
            res <- RangedData()
            for(x in 1:length(starts)) {                
                if(length(extraParams)>1) {
                    res <- rbind(res,getFeatureCounts(sites.rd[starts[x]:stops[x],], features.rd, colnam, chromSizes, widths, weightsColname, parallel=parallel, eval(as.expression(extraParams))))
                } else {
                    res <- rbind(res,getFeatureCounts(sites.rd[starts[x]:stops[x],], features.rd, colnam, chromSizes, widths, weightsColname, parallel=parallel))
                }
            }
            return(res)
        } else {
            weighted<-ifelse(is.null(weightsColname),FALSE,TRUE)
            
            # only get labels if not supplied
            if(is.null(names(widths))) {
                names(widths) <- getWindowLabel(widths)
        	}
        	
            ## use only chromosomes that are present in both sites.rd and features.rd ##
            ok.chrs <- intersect(space(sites.rd),space(features.rd))
            features.rd <- features.rd[names(features.rd) %in% ok.chrs]
            
            query <- sites.rd[,-c(1:length(colnames(sites.rd)))]
            
            if(!parallel) { registerDoSEQ() }
            
            allcounts <- foreach(x=iter(widths),.inorder=TRUE,.packages="IRanges",.export=c("query", "features.rd", "weighted", "weightsColname", "chromSizes","resizeRangedData")) %dopar% {            
                query <- resizeRangedData(query,width=x,spaceSizes=chromSizes)
                if (weighted) {
                    res <- as.data.frame(as.matrix(findOverlaps(query,features.rd,...)))
                    res$weights <- features.rd[res$subject,][[weightsColname]]
                    tapply(res$weights,res$query,sum)            
                } else {
                    ## dont use countOverlaps() since it returns overlapping ranges from other spaces/chrs if it was a factor
                    as.table(findOverlaps(query,features.rd,...))  
                }
            }
            names(allcounts)<-names(widths)
     
            for (windowName in names(widths)) {
                columnName <- paste(colnam,names(widths[ windowName ] ),sep=".")
                sites.rd[[columnName]] <- 0
                sites.rd[[columnName]][as.numeric(names(allcounts[[windowName]]))] <- as.numeric(allcounts[[windowName]])
            }
            sites.rd
        }
    }
}

#' Get counts of annotation within a defined window around each query range/position for large annotation objects spanning greater than 100 million rows. This is still in beta phase. 
#'
#' Given a query object and window size(s), the function first augments the ranges within the query by flanking starts and stops with window width. Therefore, a start of 12 and end of 14 with width 10 will yield a range of 8,17. This new range is then compared against the midpoints of subject to find any overlapping ranges and then tallied up. Note that here counting is doing using midpoint of the ranges in subject instead of start-stop boundaries. 
#'
#' @param sites.rd RangedData object to be used as the query.
#' @param features.rd RangedData obj to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves as a prefix to windows sizes!
#' @param widths a named/numeric vector of window sizes to be used for casting a net around each position. Default: \code{c(1000,10000,1000000)}
#'
#' @return a RangedData object with new annotation columns appended at the end of sites.rd. There will be a column for each width defined in widths parameter. If widths was a named vector i.e. c("100bp"=100,"1K"=1000), then the colname parameter will be pasted together with width name else default name will be generated by the function.
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{getNearestFeature}}, \code{\link{getSitesInFeature}}, \code{\link{getFeatureCounts}}.
getFeatureCountsBig <- function(sites.rd, features.rd, colnam=NULL, widths=c(1000,10000,1000000)) {
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
    
    for (windowName in names(widths)) {
        cat(".")
        columnName <- paste(colnam,names(widths[ windowName ] ),sep=".")
        sites.rd[[columnName]] <- 0
        
        query <- sites.rd[,-c(1:length(colnames(sites.rd)))]
        query <- ranges(resizeRangedData(query,width=widths[windowName],limitLess=TRUE))
        counts <- sapply(ok.chrs, function(x) {
                            features.mid<-sort((start(features.rd[[x]])+end(features.rd[[x]]))/2)                            
                            findInterval(end(query[[x]]),features.mid) - findInterval(start(query[[x]]),features.mid)
                        })             
        stopifnot(nrow(sites.rd)==sum(unlist(lapply(counts,length)))) 
        sites.rd[[columnName]][good.rows] <- unsplit(counts, space(sites.rd)[good.rows], drop = TRUE)
    }
    sites.rd
}

#' Find overlapping positions/ranges that match between the query and subject. 
#'
#' When used in genomic context, the function annotates genomic positions of interest with information like if they were in a gene or cpg island or whatever annotation that was supplied in the subject.
#'
#' @param sites.rd RangedData object to be used as the query.
#' @param features.rd RangedData obj to be used as the subject or the annotation table.
#' @param colnam column name to be added to sites.rd for the newly calculated annotation...serves a core!
#' @param asBool Flag indicating whether to return results as TRUE/FALSE or the property of an overlapping feature..namely feature name and orientation if available. Defaults to FALSE.
#' @param feature.colnam column name from features.rd to be used for retrieving the feature name. By default this is NULL assuming that features.rd has a column that includes the word 'name' somewhere in it. Not required if asBool=TRUE.
#' @param strand.colnam column name from features.rd to be used for retrieving the feature orientation. By default this is NULL assuming that features.rd has a column that includes the word 'strand' somewhere in it. If it doesn't the function will assume the supplied annotation has no orientation information '*'. Not required if asBool=TRUE.
#' @param parallel use parallel backend to perform calculation with \code{\link{foreach}}. Defaults to FALSE. Not applicable when asBool=T. If no parallel backend is registered, then a serial version of foreach is ran using \code{\link{registerDoSEQ()}}.
#' @param ... Additional parameters for \code{\link{findOverlaps}}.
#'
#' @return a RangedData object with new annotation columns appended at the end of sites.rd.
#'
#' @note Try not to use this function for >50 spaces unless you have tons fo memory. If parallel=TRUE, then be sure to have a paralle backend registered before running the function. One can use any of the following libraries compatible with \code{\link{foreach}}: doMC, doSMP, doSNOW, doMPI. For example: library(doSMP); w <- startWorkers(2); registerDoSMP(w)
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{getFeatureCounts}}, \code{\link{getNearestFeature}}.
#'
#' @examples
#' # Convert a data frame to RangedData object
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
        res <- as.data.frame(as.matrix(findOverlaps(sites.rd,features.rd,...)))
        stopifnot(!any(is.na(res)))                    

        if(is.null(feature.colnam)) {
            featureName<-getRelevantCol(colnames(features.rd),c("name","featureName"),"featureName",multiple.ok=TRUE)
            feature.colnam<-colnames(features.rd)[featureName][1]
        }
        
        skipStrandCalc<-FALSE
        if(is.null(strand.colnam)) {
            strandCol<-grep('strand',colnames(features.rd),ignore.case=TRUE,value=TRUE)[1]
            if(is.na(strandCol)) {
                skipStrandCalc<-TRUE
            } else {
                strand.colnam<-strandCol
            }
        }
                
        res$features <- unlist(values(features.rd))[feature.colnam][,1][res$subject]
        if(skipStrandCalc) {
            res$strands <- "*"
        } else {
            res$strands <- unlist(values(features.rd))[strand.colnam][,1][res$subject]
        }
        
        ## collapse rows where query returned two hits with the same featureNames due to alternative splicing or something else.
        res <- unique(res[,c("query","features","strands")]) 
        
        counts <- table(res$query)    
        res$multipleHit <- res$query %in% as.numeric(names(counts[counts>1]))
        res$featureName <- res$features
        res$strand <- res$strands
        
        if(!parallel) { registerDoSEQ() }
        
        if(any(res$multipleHit)) {
            res.unique <- unique(subset(res,!multipleHit,select=c("query","featureName","strand"),drop=TRUE))
            res <- subset(res,multipleHit,drop=TRUE)        
            
            combos <- split(res$features,res$query)
            tmp <- foreach(x=iter(combos),.inorder=TRUE) %dopar% { paste(x,collapse=",") }
            names(tmp)<-names(combos)        
            res[,"featureName"] <- unlist(tmp[as.character(res$query)])
            
            if(skipStrandCalc) {
                res[,"strand"] <- "*"
            } else {
                combos <- split(res$strands,res$query)
                tmp <- foreach(x=iter(combos),.inorder=TRUE) %dopar% { paste(x,collapse=",") }
                names(tmp)<-names(combos)        
                res[,"strand"] <- unlist(tmp[as.character(res$query)])
            }
            
            res <- unique(res[,c("query","featureName","strand")])
            res <- rbind(res,res.unique)
        } else {
            res <- unique(subset(res,!multipleHit,select=c("query","featureName","strand"),drop=TRUE))
        }
        
        ## for collapsed rows take the uniques and add them back to sites.rd        
        sites.rd[[colnam]] <- FALSE
        sites.rd[[paste(colnam,"Ort",sep="")]] <- NA
        sites.rd[[colnam]][res$query] <- as.character(res$featureName)
        sites.rd[[paste(colnam,"Ort",sep="")]][res$query] <- as.character(res$strand)
    }
    sites.rd
}

#' Annotate a RangedData object using one of annotation functions. 
#'
#' This is a wrapper function which calls one of following functions depending on annotType parameter: \code{\link{getFeatureCounts}}, \code{\link{getNearestFeature}}, code{\link{getSitesInFeature}} 
#'
#' @param annotType one of following: within, nearest, counts.
#' @param ... Additional parameters to be passed to the respective annotation function.
#'
#' @return a RangedData object with new annotation columns appended at the end of sites.rd.
#'
#' @seealso \code{\link{makeRangedData}}, \code{\link{getFeatureCounts}}, \code{\link{getNearestFeature}}, code{\link{getSitesInFeature}}.
#'
#' @examples
#' # Convert a data frame to RangedData object
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
#' doAnnotation(annotType="counts",alldata.rd,genes.rd,"NumOfGene",seqlengths(Hsapiens))
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