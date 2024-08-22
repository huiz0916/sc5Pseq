# sc5Pseq

## Script
Some useful scripts for sc5Pseq project.

### 1. `stat_base_content.py`
Description:
This script calculates and plots base frequencies in a FASTQ file and specifies bases at specific positions simultaneously.

**Usage:**
```bash
stat_base_content.py [-h] [-i INPUT_FILE] [-n NUM_BASES]
                     [-t OUTPUT_TABLE] -p OUTPUT_PLOT
                     [-s SPECIFIC_POSITIONS [SPECIFIC_POSITIONS ...]]
                     [-sp SPECIFIC_PLOT_NAME] [-a TABLE_INPUT]
```


### 2.  `star_format.py`
Description:
Parses the .final.out file and extracts the relevant information from the STAR alignment result.

**Usage:**
```bash
python star_format.py DirectoryPath OutputName
```
