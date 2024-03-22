logger::log_threshold("ERROR")

testData <- file.path(getwd(), "data")

# test location name normalization
test_that("location strings properly normalized", {
  locationsPath <- file.path(testData, "locations.rds")
  locations <- readRDS(locationsPath)
  
  for (location in names(locations)) {
    normLocation <- normalizeLocation(location)
    expect_match(locations[[location]], normLocation)
  }
})

# test parsing of GenBank files included with package
# against expected data
extData <- system.file("extdata", package="PACVr")
allFiles <- list.files(extData, recursive = TRUE)

## GenBank files
gbkTestFiles <- allFiles[grepl("(?<!-PACVr)\\.gb$", allFiles, perl = TRUE)]
isIRChecks <- c(TRUE, FALSE)

## DF files
dfData <- file.path(testData, "gbkDataDF")
dfFiles <- list.files(dfData)

## Create tests for which test data exists 
testParams <- expand.grid(gbkFile = gbkTestFiles, 
                          isIRCheck = isIRChecks,
                          stringsAsFactors = FALSE) %>%
  dplyr::mutate(dfFile = paste0(gsub("^(\\w+).+", "\\1", gbkFile), 
                                "-", isIRCheck, ".rds")) %>%
  dplyr::mutate(dfFile = dplyr::if_else(dfFile %in% dfFiles, dfFile, NA)) %>%
  dplyr::filter(!is.null(dfFile))

## Test parsing is as expected
for (index in 1:nrow(testParams)) {
  testParam <- testParams[index, ]
  dfFile <- testParam[["dfFile"]]
  gbkFile <- testParam[["gbkFile"]]
  analysisSpecs <- list(
    isIRCheck = testParam[["isIRCheck"]]
  )
  
  test_that(paste0("successful parsing of `", gbkFile, "`"), {
    gbkFileFull <- file.path(extData, gbkFile)
    dfFileFull <- file.path(dfData, dfFile)
    read.gbData <- PACVr.read.gb(gbkFileFull)
    gbkDataDF <- read.gbSeqFeatures(read.gbData, analysisSpecs)
    gbkTestDF <- readRDS(dfFileFull)
    expect_identical(gbkDataDF, gbkTestDF)
  })
}
