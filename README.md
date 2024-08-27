# sc5Pseq

## Script
Some useful scripts for sc5Pseq project.

### 1. `stat_base_content.py`
Description:  
This script calculates and plots base frequencies from a FASTQ file. It can generate an interactive HTML plot. It can also specify base positions of interest.

**Requirements**  

It could run in python2+ and python3+ enviroment, but need some libraries.  

```bash
pip3 install pandas matplotlib plotly Bio argparse gzip os
```

**Usage:**
```
usage: stat_base_content.py [-h] [-i INPUT_FILE] [-a TABLE_INPUT] -o OUTPUT_PREFIX
                            [-n NUM_BASES] [-s POS [POS ...]]
                            [-sp SPECIFIC_OUTPUT_PREFIX] [--interactive]
                            [--add_percentage]

Calculate and plot base frequencies in a FASTQ file or from an existing CSV table.

options:
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the FASTQ file, support .gz
  -a TABLE_INPUT, --table_input TABLE_INPUT
                        Path to an existing base frequencies table to use as input
                        instead of a FASTQ file.
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Output file name prefix.
  -n NUM_BASES, --num_bases NUM_BASES
                        Number of bases to analyze (default: 85).
  -s POS [POS ...], --specific_positions POS [POS ...]
                        Specific base positions or range to analyze (e.g., -s 10 15)
  -sp SPECIFIC_OUTPUT_PREFIX, --specific_output_prefix SPECIFIC_OUTPUT_PREFIX
                        Output file name prefix for specific positions plot.
  --interactive         Generate an interactive HTML plot.
  --add_percentage      Add percentage labels to the static image.
```
**Example**  
```
    module load python #(python/3.12.1)  
    - To analyze a FASTQ file and generate output files:  
        python script.py -i example.fastq -o output_prefix  
    - Specify bases at specific positions, generate dynamic interactive HTML, and add base percentage to PNG:  
        python script.py -i example.fastq -o output_prefix -s 8 16 --interactive --add_percentage -sp special_output_prefix  
    - Generate plots directly from an existing CSV table:  
        python script.py --table_input existing_frequencies.csv -o output_prefix --interactive  
```

### 2.  `star_format.py`
Description:  
Parses the .final.out file and extracts the relevant information from the STAR alignment result.

**Usage:**
```bash
python star_format.py DirectoryPath OutputName
```
