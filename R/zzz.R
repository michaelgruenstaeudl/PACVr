.onAttach = function(libname = find.package("PACViR"), pkgname = "PACViR") {
    ## CRAN packages
    #if (!require("pacman"))
    #    install.packages("pacman")
    #pacman::p_load("RCircos", "optparse", install=TRUE)
    ## Bioconductor packages
    #if (!requireNamespace("BiocManager", quietly=TRUE))
    #    install.packages("BiocManager")
    #if (!requireNamespace("genbankr", quietly=TRUE))
    #    BiocManager::install("genbankr")
}

.onLoad = function(libname = find.package("PACViR"), pkgname = "PACViR") {
    # The RCircos.Env object must be located in the global namespace.
    #assign("RCircos.Env", RCircos::RCircos.Env, .GlobalEnv)
}
