.onAttach = function(libname = find.package("PACVr"), pkgname = "PACVr") {
    # nothing
}

.onLoad = function(libname = find.package("PACVr"), pkgname = "PACVr") {
    # nothing
    require(RCircos)
}
