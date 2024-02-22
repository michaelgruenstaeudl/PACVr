logger::log_threshold("ERROR")

# test the usage files used within README.md
usageData <- system.file("extdata/README_USAGE", package="PACVr")
usageFiles <- list.files(usageData)
for (usageFile in usageFiles) {
  usageFileFull <- file.path(usageData, usageFile)
  test_that(paste0("successful run of `", usageFile, "`"), {
    returned <- suppressMessages(source(usageFileFull))
    expect_true(all(returned$value == 0))
  })
}
