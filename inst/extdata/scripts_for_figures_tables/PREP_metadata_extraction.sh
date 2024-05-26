#!/usr/bin/env bash

#contributors=c("Gregory Smith", "Nils Jenke", "Michael Gruenstaeudl")
#email="m_gruenstaeudl@fhsu.edu"
#version="2024.05.19.1300"


# this references a script that utilizes PACVr, so the package must be installed
PACVR="$(dirname "$(dirname "$(readlink -f "$0")")")/PACVr_Rscript.R"

# DEFINITIONS
SAMPLELOC=$1
ACCESSION=$2
SRA=$3

METADATA_DIR=$SAMPLELOC/${ACCESSION}.tmp
METADATA_FILE1=$METADATA_DIR/${ACCESSION}_meta_SequencingMethod.csv
METADATA_FILE2=$METADATA_DIR/${ACCESSION}_meta_AssemblyMethod.csv
METADATA_OUTF=$METADATA_DIR/${ACCESSION}_metadata.csv

# DECISION IF PACVR ALREADY RUN
echo ""
echo "Processing sample $2"
echo ""

# Test if the file $METADATA_OUTF exists
if [ -e "$METADATA_OUTF" ]; then
    echo "Metadata already extracted. Skipping PACVr ..."
	echo ""
else 
	# conduct PACVr Rscript
	echo "Running PACVr"
	echo ""
	Rscript $PACVR \
		-k $SAMPLELOC/${ACCESSION}_annotated.gb \
		-b $SAMPLELOC/${ACCESSION}_mapping_OneMoreLocations.sorted.bam \
		-t 0.5 \
		-r TRUE \
		-c TRUE \
		-o $SAMPLELOC/${ACCESSION}_CoverageDepth.pdf


	# merge metadata to one csv
	echo ""
	echo "Extracting metadata"

	esearch -db sra -query $SRA | \
		efetch -format runinfo | \
		cut -f20,7 -s -d, > $METADATA_FILE1
	sed -i 's/Model/SequencingMethod/g' $METADATA_FILE1

	esearch -db nuccore -query $ACCESSION | \
		efetch -format gb | \
		grep -E "Assembly|Sequencing" | \
		grep -v '##' | \
		sed -n '/Assembly/,/Sequencing/{/Sequencing/!p;}' | \
		awk '{$1=$1;print}' | \
		sed 's/:: /\n/g' | \
		paste -s -d ' ' | \
		sed 's/  /\n/g' | \
		awk /./ | \
		awk '{print $1}' > $METADATA_FILE2
	sed -i 's/Assembly/AssemblyMethod/g' $METADATA_FILE2
	sed -i 's/AssemblyMethod Method/AssemblyMethod/g' $METADATA_FILE2
	sed -i 's/Assembly Method/AssemblyMethod/g' $METADATA_FILE2

	if [ ! -s "$METADATA_FILE2" ]; then
		echo "Generating empty assembly metadata file ..."
		echo ""
		echo -e "AssemblyMethod\n" > $METADATA_FILE2
	fi

	find $METADATA_DIR -type f -size 0 -exec rm {} \;

	paste -d, $METADATA_DIR/*_meta_*.csv > $METADATA_OUTF
	rm $METADATA_FILE1
	rm $METADATA_FILE2

	echo ""
fi 


