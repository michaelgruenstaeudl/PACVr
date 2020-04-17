*PACVr*
=======

Plastome Assembly Coverage Visualization in R

## INSTALLATION
```
library(devtools)
install_github("michaelgruenstaeudl/PACVr")
```
Note: Detailed installation instructions can be found in the package vignette.

## USAGE
```
# In R:
library(PACVr)
gbkFile <- system.file("extdata", "NC_045072/NC_045072.gb", package="PACVr")
bamFile <- system.file("extdata", "NC_045072/NC_045072_PlastomeReadsOnly.sorted.bam", 
                       package="PACVr")
outFile <- paste(tempdir(), "/NC_045072_AssemblyCoverage_viz.pdf", sep="")
PACVr.complete(gbk.file=gbkFile, bam.file=bamFile, windowSize=250, 
               mosdepthCmd='mosdepth', logScale=FALSE, threshold=0.5,
               syntenyLineType=3, relative=TRUE, textSize=0.5,
               delete=TRUE, output=outFile)
```

## OUTPUT
![](NC_045072__all_reads.png)

## CITATION
Using PACVr in your research? Please cite it!

- Gruenstaeudl M., Jenke N. (2019). foo bar baz

```
@article {GruenstaeudlAndJenke2020,
    author = {Gruenstaeudl, M. and Jenke, N.},
    title = {PACVr: plastome assembly coverage visualization in R},
    year = {2020},
    doi = {10.1186/s12859-020-3475-0},
    URL = {https://doi.org/10.1186/s12859-020-3475-0},
    journal = {BMC Bioinformatics}
}
```

<!--
## TO DO
* Foo bar baz
* Foo bar baz
-->


## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

