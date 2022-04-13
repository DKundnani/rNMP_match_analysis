# rNMP_match_analysis
Check whether incorporated rNMPs in DNA match reference genome
Table of Contents
* [Installation](#Installation)
  * [Getting the code](#Getting the code)
  * [Creating the enviroment with required dependencies](#Creating the enviroment with required dependencies)
* [Running Mismatch Removal](#Running Mismatch Removal)
  * [Configure run](#Configure run)
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

