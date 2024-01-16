#!/usr/bin/env RScript
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.01.07.2200"

read.gbWithHandling <- function(gbkFile) {
  gbkData <- tryCatch({
    read.gb::read.gb(gbkFile, DNA=TRUE, Type="full", Source="File")
  },
    error = function(e) {
      if (conditionMessage(e) == "dim(X) must have a positive length") {
        logger::log_error(
          paste("read.gb encountered an unqualified feature;",
                "fix w/",
                "devtools::install_github(\"alephnull7/read.gb\")")
        )
      } else {
        print(e)
        logger::log_error(
          logger::skip_formatter(
            paste("read.gb encountered an unexpected error:", e)
          )
        )
      }
      return(NULL)
    })
  return(gbkData)
}
