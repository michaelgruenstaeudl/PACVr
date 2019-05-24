.onAttach = function(libname = find.package("PACViR"), pkgname = "PACViR") {
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load("RCircos", "optparse")
    # Bioconductor packages
    if (!require("BiocManager")) install.packages("BiocManager")
    source("https://bioconductor.org/biocLite.R")
    pacman::p_load("BiocParallel", "genbankr")
}

.onLoad = function(libname = find.package("PACViR"), pkgname = "PACViR") {
    pacman::p_load("RCircos", "optparse", "genbankr")
}
