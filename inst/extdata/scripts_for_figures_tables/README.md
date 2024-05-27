## Scripts for tables and figures
Scripts for generating the tables and figures of the manuscript Jenke et al. (2024)
---

### Prerequisites
- R
- R package 'PACVr' (at least version 1.1.1)
- Access to a Bash-supported terminal or terminal emulator
- [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

### Instructions
- Download the sample files of the analysis from [Zenodo](https://zenodo.org/records/4555956).
- Ensure that the Entrez Direct tools `esearch` and `efetch` are in your Bash PATH (i.e., accessible without specification of the absolute file path); if necessary, append to your PATH or set symlinks.
- With `depth_and_evenness_paper` as the current directory, execute the relevant scripts (e.g., `Rscript JenkeEtAl2024__Figure1A.R`).

### Notes

#### How to calculate the average coding/noncoding ratio for all plastid genomes under study
```
cd  folder_with_tmp_folders_from_PREP_metadata_extraction.sh

echo "" > average_coding_noncoding_length_ratio.txt
for ACC in $(ls -d *.tmp | sed 's/.tmp//g'); do
  LEN_CODING=$(cat ${ACC}.tmp/${ACC}.1_coverage.summary.genes.tsv | grep "^Unpartitioned" | awk '{print $3}')
  LEN_NONCODING=$(cat ${ACC}.tmp/${ACC}.1_coverage.summary.noncoding.tsv | grep "^Unpartitioned" | awk '{print $3}')
  echo $LEN_CODING/$LEN_NONCODING | bc -l >> average_coding_noncoding_length_ratio.txt
done
cat average_coding_noncoding_length_ratio.txt | sed '/^[[:space:]]*$/d' | paste -sd+ | bc -l
cat average_coding_noncoding_length_ratio.txt | sed '/^[[:space:]]*$/d' | wc -l
```
