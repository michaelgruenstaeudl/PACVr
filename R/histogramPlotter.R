#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.07.09.1900"

# The following R functions were taken from the R package RCircos and then modified.
# The modifications were necessary to fix several issues in the original package code.

PACVr.Get.Start.End.Locations <- function(plot.data, plot.width){
  RCircos.Cyto <- RCircos::RCircos.Get.Plot.Ideogram();

  dataChroms <- as.character(plot.data[,1]);
  chromosomes <- unique(dataChroms);
  cyto.chroms <- as.character(RCircos.Cyto$Chromosome);

  point.loc <- as.numeric(plot.data$Location)
  locations <- cbind(point.loc-plot.width, point.loc+plot.width)

  for(aChr in seq_len(length(chromosomes))){
    cyto.rows <- which(cyto.chroms==chromosomes[aChr]);
    chr.start <- min(RCircos.Cyto$StartPoint[cyto.rows]);
    chr.end   <- max(RCircos.Cyto$EndPoint[cyto.rows]);

    data.rows <- which(dataChroms==chromosomes[aChr]);
  
    start.outliers <- which(locations[data.rows, 1] < chr.start)
    which(locations[data.rows, 1] < chr.start)
    if(length(start.outliers)>0) 
      locations[data.rows[start.outliers], 1] <- chr.start;
  
    end.outliers <- which(locations[data.rows, 2] > chr.end)
    if(length(end.outliers)>0) 
      locations[data.rows[end.outliers], 2] <- chr.end;
  }
  return(locations)
}

PACVr.Histogram.Plot <- function(hist.data=NULL, data.col=4, 
                                 track.num=NULL, side=c("in", "out"), min.value=NULL, 
                                 max.value=NULL, inside.pos=NULL, outside.pos=NULL,
                                 genomic.columns=3, is.sorted=TRUE)
{
  if(is.null(hist.data)) 
    warning("Genomic data missing in input.")
    stop()
  
  boundary <- RCircos::RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                        outside.pos, FALSE);
  outerPos <- boundary[1];
  innerPos  <- boundary[2];
  
  if(is.null(genomic.columns) || genomic.columns<2 || genomic.columns>3) 
    warning("Number of columns for genomic position incorrect.")
    stop()
  if( is.null(data.col) || data.col <= genomic.columns)  
    warning(paste("Number of input columns must be > ", genomic.columns, ".", sep=""))
    stop()
  
  RCircos.Pos <- RCircos::RCircos.Get.Plot.Positions();
  RCircos.Par <- RCircos::RCircos.Get.Plot.Parameters();
  RCircos.Cyto <- RCircos::RCircos.Get.Plot.Ideogram();
  
  # Convert raw data to plot data. The raw data will be validated first during the convertion
  hist.data <- RCircos::RCircos.Get.Single.Point.Positions(hist.data,
                                                  genomic.columns);
  locations <- PACVr.Get.Start.End.Locations(hist.data, 
                                              RCircos.Par$hist.width)
  
  # Histgram colors and height
  histColors <- RCircos::RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color); 
  
  histValues <- as.numeric(hist.data[, data.col]);
  if(is.null(max.value) || is.null(min.value)) {
    max.value <- max(histValues);
    min.value <- min(histValues);
  } else {
    if(min.value > max.value) 
      warning("min.value must be greater than max.value.")
      stop()
  }
  histHeight <- RCircos::RCircos.Get.Data.Point.Height(histValues, min.value, 
                                              max.value, plot.type="points", outerPos-innerPos);
  
  # Draw histogram
  RCircos::RCircos.Track.Outline(outerPos, innerPos, RCircos.Par$sub.tracks);
  
  for(aPoint in seq_len(nrow(hist.data)))
  {
    height <- innerPos + histHeight[aPoint];
    theStart <- locations[aPoint, 1];
    theEnd <- locations[aPoint, 2];
    
  # Plot rectangle with specific height for each data point
    polygonX <- c(RCircos.Pos[theStart:theEnd,1]*height, 
                  RCircos.Pos[theEnd:theStart,1]*innerPos);
    polygonY <- c(RCircos.Pos[theStart:theEnd,2]*height, 
                  RCircos.Pos[theEnd:theStart,2]*innerPos);
    polygon(polygonX, polygonY, col=histColors[aPoint], border=NA);
  }
}
