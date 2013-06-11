#### Basic tests ####

query <- GRanges("A", IRanges(c(1, 5, 12, 20), width=1), 
                 strand=c("-","+","-","+"))
subject <- GRanges("A", IRanges(c(1,5,10,15,21), width=8:4),
                   strand=c("+", "+", "-", "-","-"))

#### test count based annotations ####
res <- findOverlaps(query, rep(subject,each=3), maxgap=10, ignore.strand=TRUE)
expect_that(countQueryHits(res), is_identical_to(c(9L, 12L, 15L, 12L)))

#### test in/out based annotations ####
res <- overlapsAny(query, subject, ignore.strand=TRUE)
expect_that(res, is_identical_to(c(TRUE, TRUE, TRUE, FALSE)))

#### test nearest type annotations ####
res <- as.data.frame(nearest(query, subject, select="all", ignore.strand=TRUE))
res <- getLowestDists(query, subject, res, "either")
rownames(res) <- NULL

## for getNearestFeature() ##
expect_that(res,
            is_identical_to(structure(list(queryHits = 1:4, 
                                           subjectHits = c(1L, 2L, 3L, 5L), 
                                           dist = c(0L, 0L, 2L, 1L)), 
                                      .Names = c("queryHits", "subjectHits", "dist"), 
                                      row.names = c(NA, -4L), class = "data.frame")))

## for get2NearestFeature() ##
res$qStrand <- as.character(strand(query))[res$queryHits]
res$u2 <- with(res, ifelse(dist<0, 
                           ifelse(qStrand=="+", subjectHits - 1, subjectHits + 1), 
                           ifelse(qStrand=="+", subjectHits - 2, subjectHits + 2)))

res$u1 <- with(res, ifelse(dist<0, subjectHits, 
                           ifelse(qStrand=="+", subjectHits - 1, subjectHits + 1)))

res$d1 <- with(res, ifelse(dist<0, 
                           ifelse(qStrand=="+", subjectHits + 1, subjectHits - 1), 
                           subjectHits))

res$d2 <- with(res, ifelse(dist<0, 
                           ifelse(qStrand=="+", subjectHits + 2, subjectHits - 2), 
                           ifelse(qStrand=="+", subjectHits + 1, subjectHits - 1)))

expect_that(res,
            is_identical_to(structure(list(queryHits = 1:4, 
                                           subjectHits = c(1L, 2L, 3L, 5L), 
                                           dist = c(0L, 0L, 2L, 1L), 
                                           qStrand = c("-", "+", "-", "+"), 
                                           u2 = c(3, 0, 5, 3), 
                                           u1 = c(2, 1, 4, 4), 
                                           d1 = c(1L, 2L, 3L, 5L), 
                                           d2 = c(0, 3, 2, 6)), 
                                      .Names = c("queryHits", "subjectHits", "dist", "qStrand", "u2", "u1", "d1", "d2"), row.names = c(NA, -4L), class = "data.frame")))

## for getNearestFeature(): 5p, 3p, and midpoint ##
subject5p <- flank(subject, width=-1) 
subject3p <- flank(subject, width=-1, start=FALSE)
ranges(subject) <- IRanges(mid(ranges(subject)), width=1)

res <- as.data.frame(nearest(query, subject5p, select="all", ignore.strand=TRUE))
res <- getLowestDists(query, subject5p, res, "5p")
expect_that(res,
            is_identical_to(structure(list(queryHits = 1:4, subjectHits = 1:4, 
                                           dist = c(0L, 0L, -3L, -1L)), 
                                      .Names = c("queryHits", "subjectHits", "dist"), 
                                      row.names = c(NA, -4L), class = "data.frame")))

res <- as.data.frame(nearest(query, subject3p, select="all", ignore.strand=TRUE))
res <- getLowestDists(query, subject3p, res, "3p")
expect_that(res,
            is_identical_to(structure(list(queryHits = 1:4, 
                                           subjectHits = c(1L, 1L, 2L, 5L), 
                                           dist = c(-7L, 3L, 1L, 1L)), 
                                      .Names = c("queryHits", "subjectHits", "dist"), 
                                      row.names = c(NA, 4L), class = "data.frame")))

res <- as.data.frame(nearest(query, subject, select="all", ignore.strand=TRUE))
res <- getLowestDists(query, subject, res, "midpoint")
expect_that(res,
            is_identical_to(structure(list(queryHits = 1:4, 
                                           subjectHits = c(1L, 1L, 3L, 5L), 
                                           dist = c(-3L, -1L, 0L, 2L)), 
                                      .Names = c("queryHits", "subjectHits", "dist"), 
                                      row.names = c(NA, 4L), class = "data.frame")))
