#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2021.04.16.2200"

# The following R functions were taken from the R package RCircos and then modified.
# The modifications were necessary to fix several issues in the original package code.


PACVr.Get.Start.End.Locations <- function(plot.data, plot.width) {
  RCircos.Cyto <- RCircos::RCircos.Get.Plot.Ideogram()
  
  
  dataChroms <- as.character(plot.data[, 1])
  
  chromosomes <- unique(dataChroms)
  
  cyto.chroms <- as.character(RCircos.Cyto$Chromosome)
  
  
  point.loc <- as.numeric(plot.data$Location)
  locations <- cbind(point.loc - plot.width, point.loc + plot.width)
  
  for (aChr in seq_len(length(chromosomes))) {
    cyto.rows <- which(cyto.chroms == chromosomes[aChr])
    
    chr.start <- min(RCircos.Cyto$StartPoint[cyto.rows])
    
    chr.end   <- max(RCircos.Cyto$EndPoint[cyto.rows])
    
    
    data.rows <- which(dataChroms == chromosomes[aChr])
    
    
    start.outliers <- which(locations[data.rows, 1] < chr.start)
    which(locations[data.rows, 1] < chr.start)
    if (length(start.outliers) > 0)
      locations[data.rows[start.outliers], 1] <- chr.start
    
    
    end.outliers <- which(locations[data.rows, 2] > chr.end)
    if (length(end.outliers) > 0)
      locations[data.rows[end.outliers], 2] <- chr.end
    
  }
  return(locations)
}


PACVr.Histogram.Plot <- function(hist.data = NULL,
                                 data.col = 4,
                                 track.num = NULL,
                                 side = c("in", "out"),
                                 min.value = NULL,
                                 max.value = NULL,
                                 inside.pos = NULL,
                                 outside.pos = NULL,
                                 genomic.columns = 3,
                                 is.sorted = TRUE)
{
  if (is.null(hist.data)) {
    warning("Genomic data missing in input.")
    stop()
  }
  boundary <-
    RCircos::RCircos.Get.Plot.Boundary(track.num, side, inside.pos,
                                       outside.pos, FALSE)
  
  outerPos <- boundary[1]
  
  innerPos  <- boundary[2]
  
  if (is.null(genomic.columns) ||
      genomic.columns < 2 || genomic.columns > 3) {
    warning("Number of columns for genomic position incorrect.")
    stop()
  }
  if (is.null(data.col) || data.col <= genomic.columns)  {
    warning(paste("Number of input columns must be > ", genomic.columns, ".", sep =
                    ""))
    stop()
  }
  RCircos.Pos <- RCircos::RCircos.Get.Plot.Positions()
  
  RCircos.Par <- RCircos::RCircos.Get.Plot.Parameters()
  
  RCircos.Cyto <- RCircos::RCircos.Get.Plot.Ideogram()
  
  # Convert raw data to plot data. The raw data will be validated first during the convertion
  hist.data <-
    RCircos::RCircos.Get.Single.Point.Positions(hist.data,
                                                genomic.columns)
  
  locations <- PACVr.Get.Start.End.Locations(hist.data,
                                             RCircos.Par$hist.width)
  # Histgram colors and height
  histColors <-
    RCircos::RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color)
  
  histValues <- as.numeric(hist.data[, data.col])
  
  if (is.null(max.value) || is.null(min.value)) {
    max.value <- max(histValues)
    
    min.value <- min(histValues)
    
  } else {
    if (min.value > max.value) {
      warning("min.value must be greater than max.value.")
      stop()
    }
  }
  histHeight <-
    RCircos::RCircos.Get.Data.Point.Height(histValues,
                                           min.value,
                                           max.value,
                                           plot.type = "points",
                                           outerPos - innerPos)
  
  # Draw histogram
  RCircos::RCircos.Track.Outline(outerPos, innerPos, RCircos.Par$sub.tracks)
  
  for (aPoint in seq_len(nrow(hist.data)))
  {
    height <- innerPos + histHeight[aPoint]
    
    theStart <- locations[aPoint, 1]
    
    theEnd <- locations[aPoint, 2]
    
    
    # Plot rectangle with specific height for each data point
    polygonX <- c(RCircos.Pos[theStart:theEnd, 1] * height,
                  RCircos.Pos[theEnd:theStart, 1] * innerPos)
    
    polygonY <- c(RCircos.Pos[theStart:theEnd, 2] * height,
                  RCircos.Pos[theEnd:theStart, 2] * innerPos)
    
    polygon(polygonX, polygonY, col = histColors[aPoint], border = NA)
    
  }
}


##### remove Chr names [UNUSED] #####
PACVr.Chromosome.Ideogram.Plot <- function(tick.interval = 0)
{
  RCircos.Draw.Chromosome.Ideogram()
  
  RCircos.Highligh.Chromosome.Ideogram()
  
  
  if (tick.interval > 0) {
    PACVr.Ideogram.Tick.Plot(tick.interval)
    
  }
}


##### tick length, text size, text orientation #####
PACVr.Ideogram.Tick.Plot <-
  function(tick.num = 10,
           track.for.ticks = 3,
           add.text.size = 0)
  {
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    
    RCircos.Pos[1:(nrow(RCircos.Pos) / 2), 3] <-
      (RCircos.Pos[1:(nrow(RCircos.Pos) / 2), 3] + 270) %% 360
    RCircos.Pos[(nrow(RCircos.Pos) / 2 + 1):nrow(RCircos.Pos), 3] <-
      (RCircos.Pos[(nrow(RCircos.Pos) / 2 + 1):nrow(RCircos.Pos), 3] + 90) %% 360
    
    endchr <- RCircos.Cyto$ChromEnd[length(RCircos.Cyto$ChromEnd)]
    tick.interval <- endchr / tick.num / 1000
    
    #   Check if there is enough space for ticks and labels. From outside
    #   of chromosome ideogram, ticks start at highlight position and take
    #   one track height, tick label takes three tracks, and chromosome
    #   names use two tracks. There will be total of 6 tracks needed.
    #   ===================================================================
    #
    track.height <- RCircos.Par$track.height
    
    tick.height  <- track.height * track.for.ticks
    
    ticks.span   <- RCircos.Par$chr.ideo.pos + tick.height * 2
    
    
    if (RCircos.Par$plot.radius < ticks.span)
    {
      stop(
        paste0(
          "There is no enough room to draw ticks.\n",
          "Please reset plot radius and redraw chromosome ideogram.\n"
        )
      )
      
    }
    
    #   Draw ticks and labels. Positions are calculated based on
    #   chromosome highlight positions
    #   ========================================================
    #
    start.pos <- RCircos.Par$highlight.pos
    
    innerPos <- RCircos.Pos[, 1:2] * start.pos
    
    outerPos <- RCircos.Pos[, 1:2] * (start.pos + track.height / 2)
    
    mid.pos  <- RCircos.Pos[, 1:2] * (start.pos + track.height / 4)
    
    
    the.interval <- tick.interval * 1000
    
    short.tick  <-
      round(the.interval / RCircos.Par$base.per.unit, digits = 0)
    
    long.tick  <-
      round(the.interval / RCircos.Par$base.per.unit, digits = 0)
    
    #short.tick <- round(the.interval/RCircos.Par$base.per.unit, digits=0);
    #long.tick  <- short.tick*2;
    
    lab.pos  <- RCircos.Pos[, 1:2] * (start.pos + tick.height / 2)
    
    chroms <- unique(RCircos.Cyto$Chromosome)
    
    for (aChr in seq_len(length(chroms)))
    {
      the.chr  <- RCircos.Cyto[RCircos.Cyto[, 1] == chroms[aChr],]
      
      chr.start <- the.chr$StartPoint[1]
      
      chr.end   <- the.chr$EndPoint[nrow(the.chr)]
      
      
      total.ticks <- tick.num
      
      for (a.tick in seq_len(total.ticks))
      {
        tick.pos <- chr.start + (a.tick - 1) * long.tick
        
        if (tick.pos < chr.end)
        {
          lines(c(innerPos[tick.pos, 1], outerPos[tick.pos, 1]),
                c(innerPos[tick.pos, 2], outerPos[tick.pos, 2]),
                col = the.chr$ChrColor[1])
          
          
          lab.text <-
            paste0(round((a.tick - 1) * tick.interval, 1), "kb")
          
          text(
            lab.pos[tick.pos, 1] ,
            lab.pos[tick.pos, 2],
            lab.text,
            cex = RCircos.Par$text.size + add.text.size,
            srt = RCircos.Pos$degree[tick.pos]
          )
          
        }
        
        tick.pos <- tick.pos + short.tick
        
        if (tick.pos < chr.end)
        {
          lines(c(innerPos[tick.pos, 1], mid.pos[tick.pos, 1]),
                c(innerPos[tick.pos, 2], mid.pos[tick.pos, 2]),
                col = the.chr$ChrColor[1])
          
        }
      }
    }
    
    #   Reset plot parameters with new chromosome name position and
    #   outside track start position. As chr.name.pos is a read-only
    #   parameter, direct work with RCircosEnvironment is needed.
    #   =======================================================
    
    old.name.pos <- RCircos.Par$chr.name.pos
    
    old.out.pos  <- RCircos.Par$track.out.start
    old.distance <- old.out.pos - old.name.pos
    
    
    RCircos.Par$chr.name.pos <- ticks.span
    
    RCircos.Par$track.out.start <-
      RCircos.Par$chr.name.pos + old.distance
    
    
    RCircosEnvironment <- NULL
    
    RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
    
    RCircosEnvironment[["RCircos.PlotPar"]] <- NULL
    
    RCircosEnvironment[["RCircos.PlotPar"]] <- RCircos.Par
    
  }


##### correction and rotate #####
PACVr.Gene.Name.Plot <- function(gene.data = NULL,
                                 name.col = NULL,
                                 track.num = NULL,
                                 side = "in",
                                 inside.pos = NULL,
                                 outside.pos = NULL,
                                 genomic.columns = 3,
                                 is.sorted = FALSE,
                                 rotate = 0,
                                 correction = 0,
                                 add.text.size = 0)
{
  if (is.null(gene.data))
    stop("Genomic data missing in RCircos.Gene.Name.Plot().\n")
  
  if (is.null(genomic.columns))
    stop("Missing number of columns for genomic position.\n")
  
  if (is.null(name.col) || name.col <= genomic.columns)
    stop("Data column must be ", genomic.columns + 1, " or bigger.\n")
  
  
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  
  textColors <-
    RCircos.Get.Plot.Colors(gene.data, RCircos.Par$text.color)
  
  RCircos.Pos[1:(nrow(RCircos.Pos) / 2), 3] <-
    (RCircos.Pos[1:(nrow(RCircos.Pos) / 2), 3] + (360 - rotate)) %% 360
  RCircos.Pos[(nrow(RCircos.Pos) / 2 + 1):nrow(RCircos.Pos), 3] <-
    (RCircos.Pos[(nrow(RCircos.Pos) / 2 + 1):nrow(RCircos.Pos), 3] + rotate) %% 360
  
  #    Convert raw data to plot data. The raw data will be validated
  #    first during the conversion
  #    =============================================================
  #
  boundary <-
    RCircos.Get.Plot.Boundary(track.num, side, inside.pos,
                              outside.pos, FALSE)
  
  gene.data <- RCircos.Get.Single.Point.Positions(gene.data,
                                                  genomic.columns)
  
  gene.data <-
    RCircos.Get.Gene.Label.Locations(gene.data,  genomic.columns,
                                     is.sorted)
  
  #    Label positions
  #    =============================================================
  #
  rightSide <- nrow(RCircos.Pos) / 2
  
  thePoints <- as.numeric(gene.data[, ncol(gene.data)])
  
  
  if (side == "in") {
    labelPos <- boundary[1] - correction
    
    textSide <- rep(4, nrow(gene.data))
    
    textSide[thePoints <= rightSide] <- 2
    
  } else {
    labelPos  <- boundary[2] - correction
    
    textSide <- rep(2, nrow(gene.data))
    
    textSide[thePoints <= rightSide] <- 4
    
  }
  
  #    Plot labels
  #    =============================================================
  #
  for (aText in seq_len(nrow(gene.data)))
  {
    geneName <- as.character(gene.data[aText, name.col])
    
    rotation <- RCircos.Pos$degree[thePoints[aText]]
    
    
    text(
      RCircos.Pos[thePoints[aText], 1] * labelPos,
      RCircos.Pos[thePoints[aText], 2] * labelPos,
      label = geneName,
      pos = textSide[aText],
      cex = RCircos.Par$text.size + add.text.size,
      srt = rotation,
      offset = 0,
      col = textColors[aText]
    )
    
  }
}


PACVr.Gene.Connector.Plot <- function(genomic.data = NULL,
                                      track.num = NULL,
                                      side = "in",
                                      inside.pos = NULL,
                                      outside.pos = NULL,
                                      genomic.columns = 3,
                                      is.sorted = FALSE)
{
  if (is.null(genomic.data))
    stop("Genomic data missing for RCircos.Gene.Connector.Plot().\n")
  
  
  boundary <-
    RCircos.Get.Plot.Boundary(track.num, side, inside.pos,
                              outside.pos, erase.area = FALSE)
  
  outerPos <- boundary[1]
  
  innerPos  <- boundary[2]
  
  
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  
  
  #    Construct Connector data from gene name data
  #    =======================================================
  #
  geneData    <- RCircos.Get.Single.Point.Positions(genomic.data,
                                                    genomic.columns)
  
  labelData   <- RCircos.Get.Gene.Label.Locations(geneData,
                                                  genomic.columns, is.sorted)
  
  connectData <-
    data.frame(labelData$Location, labelData$LabelPosition)
  
  
  #    Switch the columns of genomic location and label
  #    location for inside or outside
  #    =================================================
  #
  if (outerPos < RCircos.Par$chr.ideo.pos) {
    genomicCol <- ncol(connectData) - 1
    
    labelCol <- ncol(connectData)
    
  } else {
    genomicCol <- ncol(connectData)
    
    labelCol <- ncol(connectData) - 1
    
  }
  
  #    Heights for the two vertical lines of connectors and
  #    the horizontal line range
  #    ====================================================
  #
  vHeight <- round((outerPos - innerPos) / 10, digits = 4)
  
  hRange <- outerPos - innerPos - 2 * vHeight
  
  
  topLoc <- outerPos - vHeight
  
  botLoc <- innerPos + vHeight
  
  
  #    Connector colors
  #    ===============================================
  #
  lineColors <-
    RCircos.Get.Plot.Colors(labelData, RCircos.Par$text.color)
  
  
  #    Plot Connectors
  #    ===============================================
  #
  chroms <- unique(connectData[, 1])
  
  for (aChr in seq_along(chroms))
  {
    chrRows <- which(connectData[, 1] == chroms[aChr])
    
    total <- length(chrRows)
    
    
    for (aPoint in seq_len(total))
    {
      p1 <- connectData[chrRows[aPoint], genomicCol]
      
      p2 <- connectData[chrRows[aPoint], labelCol]
      
      
      #    draw top vertical line
      #    ======================================
      #
      lines(
        c(RCircos.Pos[p1, 1] * outerPos, RCircos.Pos[p1,  1] * topLoc),
        c(RCircos.Pos[p1, 2] * outerPos, RCircos.Pos[p1, 2] * topLoc),
        col = lineColors[chrRows[aPoint]],
        lwd = 0.1
      )
      
      
      #    draw bottom vertical line
      #    ======================================
      #
      lines(
        c(RCircos.Pos[p2, 1] * botLoc, RCircos.Pos[p2, 1] * innerPos),
        c(RCircos.Pos[p2, 2] * botLoc, RCircos.Pos[p2, 2] * innerPos),
        col = lineColors[chrRows[aPoint]],
        lwd = 0.1
      )
      
      
      #    draw horizontal line
      #    ======================================
      #
      lines(
        c(RCircos.Pos[p1,  1] * topLoc, RCircos.Pos[p2, 1] * botLoc),
        c(RCircos.Pos[p1, 2] * topLoc, RCircos.Pos[p2, 2] * botLoc),
        col = lineColors[chrRows[aPoint]],
        lwd = 0.1
      )
      
    }
  }
}

PACVr.Line.Plot <-
  function (line.data = NULL,
            data.col = 4,
            track.num = NULL,
            side = c("in", "out"),
            min.value = NULL,
            max.value = NULL,
            inside.pos = NULL,
            outside.pos = NULL,
            genomic.columns = 3,
            is.sorted = TRUE)
  {
    if (is.null(line.data))
      stop("Genomic data missing in RCircos.Line.Plot().\n")
    boundary <-
      RCircos.Get.Plot.Boundary(track.num, side, inside.pos,
                                outside.pos, FALSE)
    outerPos <- boundary[1]
    innerPos <- boundary[2]
    if (is.null(genomic.columns))
      stop("Missing number of columns for genomic position.\n")
    if (is.null(data.col) || data.col <= genomic.columns)
      stop("Line data column must be ", genomic.columns +
             1, " or bigger.\n")
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    line.data <- RCircos.Get.Single.Point.Positions(line.data,
                                                    genomic.columns)
    pointValues <- as.numeric(line.data[, data.col])
    if (is.null(min.value) || is.null(max.value)) {
      min.value <- min(pointValues)
      max.value <- max(pointValues)
    }
    else {
      if (min.value > max.value)
        stop("min.value > max.value.")
    }
    pointHeight <- RCircos.Get.Data.Point.Height(pointValues,
                                                 min.value,
                                                 max.value,
                                                 plot.type = "points",
                                                 outerPos -
                                                   innerPos)
    pointHeight <- pointHeight + innerPos
    line.colors <-
      RCircos.Get.Plot.Colors(line.data, RCircos.Par$line.color)
    #RCircos.Track.Outline(outerPos, innerPos, RCircos.Par$sub.tracks)
    for (aPoint in seq_len((nrow(line.data) - 1))) {
      if (line.data[aPoint, 1] != line.data[aPoint + 1, 1]) {
        next
      }
      point.one <- line.data[aPoint, ncol(line.data)]
      point.two <- line.data[aPoint + 1, ncol(line.data)]
      lines(
        c(RCircos.Pos[point.one, 1] * pointHeight[aPoint],
          RCircos.Pos[point.two, 1] * pointHeight[aPoint +
                                                    1]),
        c(RCircos.Pos[point.one, 2] * pointHeight[aPoint],
          RCircos.Pos[point.two, 2] * pointHeight[aPoint +
                                                    1]),
        col = line.colors[aPoint]
      )
    }
  }

PACVr.Chromosome.Ideogram.Plot <- function(tick.interval = 0)
{
  RCircos.Draw.Chromosome.Ideogram(ideo.pos = 1.125, ideo.width = 0.05)
  
  RCircos.Highligh.Chromosome.Ideogram(highlight.pos = 1.2)
  
  
  if (tick.interval > 0) {
    RCircos.Ideogram.Tick.Plot(tick.interval)
    
  }
  
  RCircos.Label.Chromosome.Names()
  
}


PACVr.Reset.Plot.Parameters <- function (new.params = NULL)
{
  if (is.null(new.params))
    stop("Missing function argument.\n")
  
  old.params <- RCircos.Get.Plot.Parameters()
  
  
  #   1.  If parameters related to total number of data tracks
  #       need reset, use reset core components instead.
  #   ==========================================================
  #if( new.params$radius.len != old.params$radius.len ||
  #    new.params$plot.radius != old.params$plot.radius ||
  #    new.params$chr.ideo.pos != old.params$chr.ideo.pos ||
  #    new.params$tracks.inside != old.params$tracks.inside ||
  #    new.params$tracks.outside != old.params$tracks.outside )
  #{ stop("Please use RCircos.Set.Core.Components() instead.\n") }
  
  #   2.  If parameters related to chromosome ideogram plot
  #       need reset, use customized plot methods instead.
  #   ==========================================================
  #if( new.params$chr.ideo.pos != old.params$chr.ideo.pos ||
  #    new.params$highlight.pos != old.params$highlight.pos ||
  #    new.params$chr.name.pos  != old.params$chr.name.pos )
  #{ stop("Please use customized ideogram plot methods instead.\n")}
  
  #   3.  Validate the parameter values in case of multiple
  #       parameters were reset in new.params such as nemeric
  #       values and color values, and ideogram layout values.
  #   ==========================================================
  RCircos.Validate.Plot.Parameters(new.params)
  
  
  #   4.  Parameters related to ideogram width change.
  #       Note: chr.ideo.pos is a read-only parameter
  #   ========================================================
  if (new.params$chrom.width != old.params$chrom.width)
  {
    differ <- new.params$chrom.width - old.params$chrom.width
    
    new.params$highlight.pos <- old.params$highlight.pos + differ
    
    
    new.name.pos <- old.params$chr.name.pos + differ
    
    if (new.params$chr.name.pos < new.name.pos)
      new.params$chr.name.pos  <- new.name.pos
    
  }
  
  #   5.  In case user modified track.in.start and track.out.start
  #   ===========================================================
  #if(new.params$track.in.start >= new.params$chr.ideo.pos)
  #  new.params$track.in.start <- new.params$chr.ideo.pos - 0.05;
  #
  #  new.name.end <- new.params$chr.name.pos + 0.3;
  #  if(new.params$track.out.start < new.name.end)
  #    new.params$track.out.start <- new.name.end;
  
  #   6.  Parameters related to data track layout. If reset, total
  #       tracks will be different. Just validate and give a prompt
  #   ============================================================
  #if(new.params$track.padding !=  old.params$track.padding ||
  #   new.params$track.height != old.params$track.height )
  #{
  #  message(paste0("Track height and/or track padding have been ",
  #                 "reset\n. Actual total data track plotted may differ.\n"));
  #}
  
  #   7.  Adjust chromosome padding parameter with default constant
  #       if base.per.unit was rest but chrom.padding was unchanged
  #       or the new chrom.paddings is too big
  #   =============================================================
  if (old.params$base.per.unit != new.params$base.per.unit &&
      old.params$chrom.paddings == new.params$chrom.paddings)
  {
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    
    band.len <- RCircos.Cyto$ChromEnd - RCircos.Cyto$ChromStart
    
    genome.len <- sum(as.numeric(band.len))
    
    padding.const <- RCircos.Get.Padding.Constant()
    
    
    total.units <- genome.len / new.params$base.per.unit
    
    new.padding <- round(padding.const * total.units, digits = 0)
    
    
    if (new.padding != new.params$base.per.unit) {
      message(
        paste(
          "\nNote: chrom.padding",
          new.params$chrom.paddings,
          " was reset to",
          new.padding,
          "\n"
        )
      )
      
      new.params$chrom.paddings <- new.padding
      
    }
  }
  
  #   Save new parameters to RCircos Environment
  #   =====================================================
  RCircosEnvironment <- NULL
  
  RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
  
  RCircosEnvironment[["RCircos.PlotPar"]] <- NULL
  
  RCircosEnvironment[["RCircos.PlotPar"]] <- new.params
  
  
  #   Ideogram/band positions are binded to base.per.unit and
  #   chromosome padding so have to be reset if base.per.unit
  #   and/or chrom.paddings are reset.
  #   ======================================================
  if (old.params$base.per.unit != new.params$base.per.unit ||
      old.params$chrom.paddings != new.params$chrom.paddings)
  {
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    
    RCircos.Cyto <- RCircos.Cyto[, 1:5]
    
    
    RCircosEnvironment[["RCircos.Cytoband"]] <- NULL
    RCircos.Set.Cytoband.Data(RCircos.Cyto)
    
    
    RCircosEnvironment[["RCircos.Base.Position"]] <- NULL
    RCircos.Set.Base.Plot.Positions()
    
  }
}
