.onAttach = function(libname = find.package("PACVr"), pkgname = "PACVr") {
    # nothing
}

.onLoad = function(libname = find.package("PACVr"), pkgname = "PACVr") {
    # nothing
    library(dplyr)
    library(logger)
    library(RCircos)
}
