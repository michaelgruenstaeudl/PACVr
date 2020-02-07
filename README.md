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
gbkFile <- system.file("extdata", "MH161174/MH161174.gb", package="PACVr")
bamFile <- system.file("extdata", "MH161174/MH161174_PlastomeReadsOnly.sorted.bam", 
                       package="PACVr")
outFile <- paste(tempdir(), "/MH161174_AssemblyCoverage_viz.pdf", sep="")
PACVr.complete(gbk.file=gbkFile, bam.file=bamFile, windowSize=250, 
               mosdepthCmd='mosdepth', logScale=FALSE, threshold=0.5,
               syntenyLineType=3, relative=TRUE, textSize=0.5,
               delete=TRUE, output=outFile)
```

## OUTPUT
![](MH161174_AssemblyCoverage_viz.png)

<!--
## CITATION
Using PACVr in your research? Please cite it!

- Gruenstaeudl M., Jenke N. (2019). foo bar baz

```
@article {Gruenstaeudl435644,
    author = {Gruenstaeudl, Michael and Hartmaring, Yannick},
    title = {EMBL2checklists: A Python package to facilitate the user-friendly submission of plant DNA barcoding sequences to ENA},
    elocation-id = {435644},
    year = {2018},
    doi = {10.1101/435644},
    URL = {https://www.biorxiv.org/content/early/2018/10/05/435644},
    journal = {bioRxiv}
}
```
-->


<!--
## TO DO
* Foo bar baz
* Foo bar baz
-->


## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

