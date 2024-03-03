*PACVr*
=======
<!-- badges: start -->
[![prepare-release](../../actions/workflows/prepare-release.yaml/badge.svg)](../../actions/workflows/prepare-release.yaml)
<!-- badges: end -->

Plastome Assembly Coverage Visualization in R

## INSTALLATION
```
library(devtools)
install_github("michaelgruenstaeudl/PACVr")
```
Note: Detailed installation instructions can be found in the package vignette.

## USAGE
### Plastid Genome
```
# In R:
library(PACVr)
gbkFile <- system.file("extdata", "NC_045072/NC_045072.gb", package="PACVr")
bamFile <- system.file("extdata", "NC_045072/NC_045072_subsampled.bam", 
                       package="PACVr")

outFile <- paste(tempdir(), "/NC_045072_CoverageViz.pdf", sep="")
#outFile <- "../NC_045072_CoverageViz.pdf"  # on R-Studio for Windows
#outFile <- "~/NC_045072_CoverageViz.pdf"   # on R-Studio for Linux

## ONLY COVERAGE VALUES, NO REGION INDICATORS ##
PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE, 
               threshold=0.5, relative=TRUE, textSize=0.5, 
               output=outFile)

## COVERAGE VALUES PLUS REGION INDICATORS ##
PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE, 
               threshold=0.5, relative=TRUE, textSize=0.5, 
               IRCheck=0, output=outFile)

## COVERAGE VALUES PLUS REGION INDICATORS PLUS IR SYNTENY LINES ##
PACVr.complete(gbkFile, bamFile, windowSize=250, logScale=FALSE, 
               threshold=0.5, relative=TRUE, textSize=0.5, 
               IRCheck=1, output=outFile)
```

### Mitochondrial Genome
```
# In R:
library(PACVr)
gbkFile <- system.file("extdata", "MG936619/MG936619.gb", package="PACVr")
bamFile <- system.file("extdata", "MG936619/MG936619_subsampled.bam", 
                       package="PACVr")

outFile <- paste(tempdir(), "/MG936619_CoverageViz.pdf", sep="")
#outFile <- "../MG936619_CoverageViz.pdf"  # on R-Studio for Windows
#outFile <- "~/MG936619_CoverageViz.pdf"   # on R-Studio for Linux

## ONLY COVERAGE VALUES, NO REGION INDICATORS ##
PACVr.complete(gbkFile, bamFile, windowSize=50, logScale=FALSE, 
               threshold=0.5, relative=TRUE, textSize=0.5, 
               output=outFile)

```

## OUTPUT
![](NC_045072_CoverageViz.png)

## FULLY AUTOMATED
```
SMPLNME="Cb01A_IOGA"
INFASTA=${SMPLNME}.fasta
READSR1=Cb01A_PlastomeReadsOnly_R1.fastq
READSR2=Cb01A_PlastomeReadsOnly_R2.fastq

mkdir -p db
bowtie2-build $INFASTA db/${SMPLNME}
bowtie2 -x db/${SMPLNME} -1 $READSR1 -2 $READSR2 -S ${SMPLNME}_mapping.sam
samtools view -Sb -F 0x04 ${SMPLNME}_mapping.sam > ${SMPLNME}_mapping_OneMoreLocations.bam
samtools sort ${SMPLNME}_mapping_OneMoreLocations.bam > ${SMPLNME}_mapping_OneMoreLocations.sorted.bam
rm $(ls *.?am | grep -v sorted)
samtools index ${SMPLNME}_mapping_OneMoreLocations.sorted.bam
# samtools view -s 0.5 -b ${SMPLNME}_mapping_OneMoreLocations.sorted.bam | samtools view -bS > ${SMPLNME}_subsampled.bam  # subsampling, if desired
Rscript PACVr/inst/extdata/PACVr_Rscript.R -k ${SMPLNME}.gb -b ${SMPLNME}_mapping_OneMoreLocations.sorted.bam -t 0.5 -r TRUE -d FALSE -o ${SMPLNME}_CoverageDepth.pdf
```

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

<!--
## PRE-FORMATTING INPUT
Due to the internal usage of R package [genbankr](https://bioconductor.org/packages/release/bioc/html/genbankr.html), any GenBank flatfile must conform to the following specifications: 
- Flatfile must include a _source_ feature at start of feature table
- All _exon_ features (plus their qualifier lines) must be removed: `sed -i -e '/    exon/,+2d' input.gb`
- All redundant _complement_ specifications must be removed: `sed -i -z 's/),\s*complement(/,/g' input.gb`
- All duplicate lines, if any, must be removed: `sed -i '$!N; /^\(.*\)\n\1$/!P; D' input.gb`
-->

## CHANGELOG
See [`CHANGELOG.md`](CHANGELOG.md) for a list of recent changes to the software.

