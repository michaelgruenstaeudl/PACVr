logger::log_threshold("ERROR")

extData <- system.file("extdata", package="PACVr")
allFiles <- list.files(extData, recursive = TRUE)
gbkTestFiles <- allFiles[grepl("(?<!-PACVr)\\.gb$", allFiles, perl = TRUE)]

for (gbkFile in gbkTestFiles) {
  test_that("successful parsing of `{gbkFile}`", {
    gbkFileFull <- file.path(extData, gbkFile)
    gbkData <- suppressMessages(PACVr.read.gb(gbkFileFull))
    gbkDataDF <- read.gb2DF(gbkData, FALSE)
    expect_false(is.null(gbkDataDF))
  })
}
