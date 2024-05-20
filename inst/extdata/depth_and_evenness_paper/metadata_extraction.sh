#!/bin/bash
#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"


# this references a script that utilizes PACVr, so the package must be installed
PACVR="$(dirname "$(dirname "$(readlink -f "$0")")")/PACVr_Rscript.R"

# this corresponds to the installation directory of Entrez
NCBIENT=/PATH/TO/NCBI-ENTREZ-DIRECT

# the sample files are not included here, but can be accessed from https://zenodo.org/records/4555956
SAMPLELOC=/PATH/TO/SAMPLES

# DEFINITIONS
ACCESSION=$1
SRA=$2

# conduct PACVr Rscript
Rscript $PACVR -k $SAMPLELOC/${ACCESSION}_annotated.gb -b $SAMPLELOC/${ACCESSION}_mapping_OneMoreLocations.sorted.bam -t 0.5 -r TRUE -c TRUE -o $SAMPLELOC/${ACCESSION}_CoverageDepth.pdf

# merge metadata to one csv
${NCBIENT}esearch -db sra -query $SRA | efetch -format runinfo | cut -f20,7 -s -d, > $SAMPLELOC/${ACCESSION}.tmp/${ACCESSION}_SRA_meta.csv
${NCBIENT}esearch -db nuccore -query $ACCESSION | efetch -format gb | sed -n '/Assembly Method/,/Sequencing/{/Sequencing/!p;}'  | awk '{$1=$1;print}' | sed 's/:: /\n/g' | paste -s -d ' ' | sed 's/  /\n/g' | awk /./ > $SAMPLELOC/${ACCESSION}.tmp/${ACCESSION}_assembly_tech.csv
find $SAMPLELOC/${ACCESSION}.tmp/ -type f -size 0 -exec rm {} \;
paste -d, $SAMPLELOC/${ACCESSION}.tmp/*.csv > $SAMPLELOC/${ACCESSION}.tmp/${ACCESSION}_metadata.csv
