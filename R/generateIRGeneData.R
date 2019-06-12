#!/usr/bin/R
#contributors = c("Michael Gruenstaeudl","Nils Jenke")
#email = "m.gruenstaeudl@fu-berlin.de", "nilsj24@zedat.fu-berlin.de"
#version = "2019.06.12.1530"

GenerateIRGeneData <- function(sample_genes) {
  # Function to generate link data for RCircos.Link.Line.Plot
  # ARGS:
  #   sample_genes: data.frame of sample data
  # RETURNS:
  #   data.frame with paired regions
  if (ncol(sample_genes) < 3) {
    stop("Dataframe must have at least 3 columns (containing chromosome name, chromosome start and chromosome end).\n")
  }
  ira <- sample_genes[sample_genes[1] == "IRa",]
  ira <- ira[nrow(ira):1,]
  ira <- ira[1:4]
  ira <- subset(ira,!grepl("IRb|IRa", ira[,"gene"],ignore.case = FALSE))
  irb <- sample_genes[sample_genes[1] == "IRb",]
  irb <- irb[1:4]
  irb <- subset(irb,!grepl("IRb|IRa", irb[,"gene"],ignore.case = FALSE))
  ira <- ira[order(ira[,2],decreasing = TRUE),]
  irb <- irb[order(irb[,3],decreasing = FALSE),]
  if (length(irb[ ,4]) != length(ira[ ,4])) {
    stop("Repeat regions need to be identical in length.")
  }
  if (!identical(ira[,4],irb[,4])) {
    stop("Repeat regions need to be identical in annotations.")
  }
  ira <- as.integer(rowMeans(ira[,2:3]))
  irb <- as.integer(rowMeans(irb[,2:3]))
  link_data <- data.frame(Chromosome = rep("IRa",length(ira)), 
                          chromStart = (ira),
                          chromEnd   = (ira),
                          Chromosome.01 = rep("IRb",length(irb)), 
                          chromStart.01 = (irb),
                          chromEnd.01   = (irb))
  return(link_data)
}
