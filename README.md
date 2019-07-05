*PACVr*
=======

Plastome Assembly Coverage Visualization in R

## INSTALLATION
Coming soon.

<!---
```
# In R:
install.packages("PACVr")
```
--->

## USAGE
```
# In R:
library(PACVr)

gbkFile <- system.file("extdata", "MH161174/MH161174.gb", package="PACVr")
bamFile <- system.file("extdata", "MH161174/MH161174_PlastomeReadsOnly.sorted.bam", 
                       package="PACVr")
outFile <- paste(getwd(), "/MH161174_AssemblyCoverage_viz.pdf", sep="")

PACVr.complete(gbk.file=gbkFile, bam.file=bamFile, windowSize=250, 
                mosdepthCmd='mosdepth', threshold=15, delete=TRUE, output=outFile)
```

## OUTPUT
![](vignettes/MH161174_AssemblyCoverage_viz.pdf)

<!---
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

<!---
## DEVELOPMENT
1. Currently nothing.
2. Currently nothing.
--->

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

