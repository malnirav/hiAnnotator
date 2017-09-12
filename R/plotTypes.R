#' Plot distance distribution to a feature boundary.
#'
#' Given a dataframe of samples and distance based annotation, the function
#' calculates the distribution of data in or around the given annotation.
#' From genomic point of view, the function can be used to identify
#' distribution of data  around genomic features like gene TSS, CpG island, etc.
#'
#' @param dat a dataframe/GRanges with required columns to make the plot.
#' @param grouping name of the column grouping the data or denoting the samples
#' @param annotCol name of the column holding the distance to feature data.
#' This can also be boolean data in which case plot will be in/out of feature.
#' @param breaks intervals by which to break up the distance data. Default is
#' seq(-1e5,1e5,5e3). Not required if `annotCol` is of type boolean.
#' @param discreteBins whether to plot continuous variable supplied in annotCol
#' as a discrete axis. This conserves plotting area, thus default is TRUE.
#' @param geom plot distribution using bars or lines? Default is 'bar'. One can
#' use 'line' as well when there are many groups.
#' @param stacked make a stacked plot? Only applies when geom is 'bar'.
#' Default is FALSE.
#' @param typeRatio whether to plot data as ratio of experimental to controls.
#' Default is FALSE. Enabling this requires a column in 'dat' called "type" with
#' two values "expr" for experimental and "ctrl" for control. This column
#' subdivides data within each group. Enabling this transforms the
#' data into plotting distribution of ratios of experimental/controls around
#' feature of interest.
#' @param printPlotData return summarized plot data? Default is FALSE.
#'
#' @return ggplot2 plot and/or table of summarized plot data.
#'
#' @seealso  \code{\link{makeGRanges}}, \code{\link{getNearestFeature}},
#' \code{\link{getSitesInFeature}}, \code{\link{getFeatureCounts}}
#'
#' @export
#'
#' @examples
#' # Convert a dataframe to GRanges object
#' data(sites)
#' data(sites.ctrl)
#' sites$type <- "expr"
#' sites <- rbind(sites,sites.ctrl)
#' alldata.rd <- makeGRanges(sites,soloStart=TRUE)
#'
#' data(genes)
#' genes.rd <- makeGRanges(genes)
#'
#' res <- doAnnotation(annotType="within", alldata.rd, genes.rd, "InGene",
#' asBool=TRUE)
#' plotdisFeature(res, "virus", "InGene")
#' plotdisFeature(res, "virus", "InGene", typeRatio=TRUE)
#' \dontrun{
#' res <- doAnnotation(annotType="nearest", res, genes.rd, "NearestGene",
#' side='5p')
#' plotdisFeature(res, "virus", "X5pNearestGeneDist")
#' plotdisFeature(res, "virus", "X5pNearestGeneDist", typeRatio=TRUE)
#' }
plotdisFeature <- function(dat = NULL, grouping = NULL, annotCol = NULL,
                           breaks = NULL, discreteBins = TRUE, geom = 'bar',
                           stacked = FALSE, typeRatio = FALSE,
                           printPlotData = FALSE) {

  ## this is to avoid "no visible binding for global variable" in R CMD check
  Distance <- type <- Percent <- Ratio <- DistToFeature <- NULL

  if(is.null(dat)) {
    stop("No data supplied")
  } else {
    if(is(dat,"GRanges")) {
      dat <- as.data.frame(dat)
    }
  }

  if(is.null(grouping)) {
    stop("No grouping parameter defined")
  } else {
    if(!grouping %in% names(dat)) {
      stop("Variable '", grouping, "' not found in supplied data!")
    }
  }
  dat$grouping <- dat[,grouping]

  if(is.null(annotCol)) {
    stop("No annotCol parameter defined")
  } else {
    if(!annotCol %in% names(dat)) {
      stop("Variable '", annotCol, "' not found in supplied data!")
    }
  }
  dat$annotCol <- dat[,annotCol]

  if(is.null(breaks)) {
    if(!all(is.logical(dat[,annotCol]))) {
      breaks <- seq(-1e5, 1e5, 5e3)
    }
  }

  if(!"type" %in% names(dat)) {
    dat$type <- ""
  }

  if(typeRatio) {
    test <- setdiff(c("expr", "ctrl"), tolower(unique(dat$type)))
    if(length(test) > 0) {
      stop("Value(s) ", paste(test, collapse = " or "),
           " not found in 'type' variable!")
    }
  }

  isBool <- FALSE
  counts.df <- dat %>% count(type, grouping)
  names(counts.df)[3] <-  "Total"
  if(all(is.logical(dat[, annotCol]))) {
    message("performing boolean summary")
    isBool <- TRUE
    plot.frame <- dat %>% count(type, annotCol, grouping) %>% rename(freq = n)
  } else {
    dat$Distance <- cut(dat$annotCol, breaks = breaks, include.lowest = TRUE,
                        dig.lab = 5)
    plot.frame <- dat %>% count(type, grouping, Distance) %>% rename(freq = n)
    plot.frame <- subset(plot.frame, !is.na(Distance))
    plot.frame$DistToFeature <- as.numeric(sub(".+,(.+)]", "\\1",
                                               plot.frame$Distance))
  }

  plot.frame <- merge(plot.frame, counts.df)
  plot.frame$Percent <- with(plot.frame, freq / Total)

  if(isBool & !discreteBins) {
    message("Data is found to be discrete not continuous.")
    discreteBins <- TRUE
  }

  if(typeRatio) {
    if(isBool) {
      plot.frame <- plot.frame %>% arrange(grouping, type, annotCol)

      ratiosCalc <- plot.frame %>% group_by(grouping, annotCol) %>%
        mutate(Ratio = Percent[tolower(type) == "expr"] /
                 Percent[tolower(type) != "expr"]) %>% ungroup %>%
        dplyr::filter(tolower(type) == "expr")

      p <- ggplot(data = ratiosCalc, aes(x = annotCol, y = Ratio)) +
        scale_x_discrete(annotCol, expand = c(0, 0))
    } else {
      plot.frame <- plot.frame %>% arrange(grouping, type, DistToFeature)

      ratiosCalc <- plot.frame %>%
        group_by(grouping, Distance, DistToFeature) %>%
        mutate(Ratio = last(Percent) / first(Percent)) %>% ungroup %>%
        dplyr::filter(tolower(type) == "expr")

      if(discreteBins) {
        p <- ggplot(data = ratiosCalc, aes(x = Distance, y = Ratio)) +
          scale_x_discrete(annotCol, expand = c(0, 0))
      } else {
        p <- ggplot(data = ratiosCalc, aes(x = DistToFeature, y = Ratio)) +
          scale_x_continuous(annotCol, expand = c(0, 0))
      }
    }

    if(geom == 'bar') {
      p <- p + geom_bar(stat = "identity", aes(fill = grouping),
                        position = ifelse(stacked, "stack", "dodge"))
    } else {
      p <- p + geom_line(aes(colour = grouping, group = grouping))
    }

    p <- p + scale_y_continuous("Ratio Expr/Ctrls", expand = c(0, 0)) +
      geom_hline(yintercept = 1)
  } else {
    if(isBool) {
      p <- ggplot(data = plot.frame, aes(x = annotCol, y = Percent)) +
        scale_x_discrete(annotCol, expand = c(0, 0))
    } else {
      if (discreteBins) {
        p <- ggplot(data = plot.frame, aes(x = Distance, y = Percent)) +
          scale_x_discrete(annotCol, expand = c(0, 0))
      } else {
        p <- ggplot(data = plot.frame, aes(x = DistToFeature, y = Percent)) +
          scale_x_continuous(annotCol, expand = c(0, 0))
      }
    }

    if(geom == 'bar') {
      p <- p + geom_bar(stat = "identity", aes(fill = grouping),
                        position = ifelse(stacked, "stack", "dodge"))
    } else {
      p <- p + geom_line(aes(colour = grouping, group = grouping))
    }

    p <- p + scale_y_continuous("Percent of Sites",
                                labels = percent, expand = c(0, 0))

    if(all(dat$type != "")) {
      p <- p + facet_wrap( ~ type, ncol = 1)
    }
  }

  p <- p + theme_bw()

  if(!isBool) {
    p <- p + theme(axis.text.x = element_text(angle = 90,
                                              hjust = 1,
                                              vjust = 0.5))
  }

  n <- length(as.character(unique(dat[, grouping])))
  allCols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
               "#A65628", "#F781BF", "#999999", "#1B9E77", "#D95F02", "#7570B3",
               "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#7FC97F",
               "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
               "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
               "#E5C494", "#B3B3B3", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
               "#FB9A99", "#E31A1C", "#FDBF6F", "#CAB2D6", "#6A3D9A", "#B15928",
               "#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598",
               "#ABDDA4", "#3288BD", "#5E4FA2")

  if(n >= 8 & n <= length(allCols)) {
    dat.colors <- as.character(allCols[1:n])
    if(geom == 'bar') {
      p <- p + scale_fill_manual(name = grouping, values = dat.colors)
    } else {
      p <- p + scale_colour_manual(name = grouping, values = dat.colors)
    }
  } else {
    if(geom == 'bar') {
      p <- p + scale_fill_discrete(name = grouping)
    } else {
      p <- p + scale_colour_discrete(name = grouping)
    }
  }

  if(printPlotData) {
    print(plot.frame)
  }

  p
}

