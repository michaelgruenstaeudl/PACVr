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
