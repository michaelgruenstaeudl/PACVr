logger::log_threshold("ERROR")

# test GenBank files included with package
extData <- system.file("extdata", package="PACVr")
allFiles <- list.files(extData, recursive = TRUE)
gbkTestFiles <- allFiles[grepl("(?<!-PACVr)\\.gb$", allFiles, perl = TRUE)]

for (gbkFile in gbkTestFiles) {
  test_that(paste0("successful parsing of `", gbkFile, "`"), {
    gbkFileFull <- file.path(extData, gbkFile)
    gbkData <- suppressMessages(PACVr.read.gb(gbkFileFull))
    gbkDataDF <- read.gb2DF(gbkData, FALSE)
    expect_false(is.null(gbkDataDF))
  })
}

# test Entrez sourced GenBank data
nucIDs <- c(
  "MG936619"
  )
for (nucID in nucIDs) {
  test_that(paste0("successful parsing of `", nucID, "`"), {
    gbkChar <- rentrez::entrez_fetch(db = "nuccore",
                                     id = nucID, 
                                     rettype = "gb")
    gbkData <- suppressMessages(PACVr.read.gb(gbkChar))
    gbkDataDF <- read.gb2DF(gbkData, FALSE)
    expect_false(is.null(gbkDataDF))
  })
}
