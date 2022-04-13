# rNMP_match_analysis
Check whether incorporated rNMPs in DNA match reference genome
Table of Contents
* [Installation](#Installation)
  * [Getting the code](#getting-the-code)
  * [Creating the enviroment with required dependencies](#Creating-the-enviroment-with-required-dependencies)
* [Running Mismatch Removal](#Running-Mismatch-Removal)
  * [Configure run](#Configure-run)
  * [Filtration](#Filtration)
* [Output](#Output)




## Installation

### Getting the code
```bash
git clone https://github.com/DKundnani/rNMP_match_analysis.git 
```

### Creating the enviroment with required dependencies
```bash
conda env create --file rNMP_match_analysis/mm_removal.yml
```

## Running Mismatch Removal
### Configure run
```bash
vim MMremoval_configure_run.sh
#Change the variables for you run as per mentions in the bash configure file
```

### Filtration
```bash
bash MMremoval_configure_run.sh

```
## Output
* From Mismatch Removal
  * Files with prefix 'MManalysis_' : chr start stop upstream+rNMP(sequencing data) upstream+rNMP(reference genome) strand (6 columns bed files)
  * Files with prefix 'poly_' : chr start stop upstream+rNMP(sequencing data) upstream+rNMP(reference genome) strand rNMP(sequencing data) rNMP(reference genome) (8 columns bed files)
  * files with prefix 'matches_' : Files with prefix 'poly_' filtered for matched column 7 & column 8
  * Files with suffix 'final.bed' : Files filtered for matched polyN upto user specified bases upstream of rNMP location
  * chrM_mismatch_percentage.txt & chrM_mismatch_percentage.txt: Difference between the percetage of rNMP bases between the bed file provided and final bed in mitochondria and nucleus

