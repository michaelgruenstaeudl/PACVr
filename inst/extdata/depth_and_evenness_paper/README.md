## Table and Figure Creation

### Depth and evenness of sequencing coverage are associated with assembly quality, genome structure, and choice of sequencing platform in archived plastid genomes

---

### Instructions

Detailed below are the methods that were used to produce many of the tables and figures 
used in the paper listed above. Reproducing these results is mostly automated, but
it does require some setup and minor changes to a couple of files. These actions are necessary
to reflect unknowns about a specific system's configuration, especially the location of
the necessary sample `.gb` and `.bam` files.

1. At minimum, have R and the `PACVr` package installed on your system, 
including all of required dependencies. Instructions for installing both R
and RStudio are located at [Posit](https://posit.co/download/rstudio-desktop/).
2. Gain access to a Bash-supported terminal emulator on your system. 
On a Unix-Like system such as Linux or macOS this should already be on your system. 
For Windows, [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) and 
[Git BASH](https://gitforwindows.org/) are recommended options.
Alternatively, you can utilize a terminal emulator included in an IDE such as RStudio.
3. Install Entrez on your system. If not already familiarized with the tool, 
[NIH Tools](https://www.ncbi.nlm.nih.gov/home/tools/) is an excellent reference.
4. Download the sample files used in the analysis. The current version of the 
archive can be found at [Zenodo](https://zenodo.org/records/4555956).
5. Within `metadata_extraction.sh`, update `NCBIENT` to the installation path of Entrez,
and update `SAMPLELOC` to the location of the uncompressed sample files.
6. Within `coverage_data_assembly.R`, in the first section `### INIT NEEDED DIRECTORIES ###`,
update `w_dir` to the location of the uncompressed sample files.
7. Using a terminal emulator with `depth_and_evenness_paper` as the current directory, execute `Rscript main.R`.
