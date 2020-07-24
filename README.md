*PACVr*
=======

Plastome Assembly Coverage Visualization in R

## INSTALLATION
```
library(devtools)
install_github("michaelgruenstaeudl/PACVr")
```
Note: Detailed installation instructions can be found in the package vignette.

## PRE-FORMATTING INPUT
Due to the internal usage of R package [genbankr](https://bioconductor.org/packages/release/bioc/html/genbankr.html), any GenBank flatfile must conform to the following specifications: 
- Flatfile must include a _source_ feature at start of feature table
- All exon features (plus their qualifier lines) must be removed: `sed -i -e '/    exon/,+2d' input.gb`
- All duplicate lines must be removed: `sed -i '$!N; /^\(.*\)\n\1$/!P; D' input.gb`
- All redundant complement specifications must be removed. For example: `join(complement(4316..4354),complement(1752..1790))  -->  join(complement(4316..4354,1752..1790))`)


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

```
@article {GruenstaeudlAndJenke2020,
    author = {Gruenstaeudl, M. and Jenke, N.},
    title = {PACVr: plastome assembly coverage visualization in R},
    year = {2020},
    doi = {10.1186/s12859-020-3475-0},
    journal = {BMC Bioinformatics},
    volume = {21},
    pages = {207}
}
```

<!--
## TO DO
* Foo bar baz
* Foo bar baz
-->


## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

