# sc5Pseq

## Script
Some useful scripts for sc5Pseq project.

### 1. `stat_base_content.py`
Description:  
This script calculates and plots base frequencies in a FASTQ file. It can also specify base positions of interest

**Usage:**
```
stat_base_content.py [-h] [-i INPUT_FILE] [-n NUM_BASES]
                            [-t OUTPUT_TABLE] -p OUTPUT_PLOT
                            [-s POS [POS ...]] [-sp SPECIFIC_PLOT_NAME]
                            [-a TABLE_INPUT]

Calculate and plot base frequencies in a FASTQ file or from an existing
frequency table.

optional arguments:
  -h, --help   show this help message and exit
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the FASTQ file, support .gz
  -n NUM_BASES, --num_bases NUM_BASES
                        Number of bases to analyze (default: 85).
  -t OUTPUT_TABLE, --output_table OUTPUT_TABLE
                        Output table file path.
  -p OUTPUT_PLOT, --output_plot OUTPUT_PLOT
                        Output plot file path.
  -s POS [POS ...], --specific_positions POS [POS ...]
                        Specific base positions to analyze (e.g., -s 10 15)
  -sp SPECIFIC_PLOT_NAME, --specific_plot_name SPECIFIC_PLOT_NAME
                        Custom name for the specific positions plot output.
  -a TABLE_INPUT, --table_input TABLE_INPUT
                        Path to an existing base frequencies table to use as
                        input instead of a FASTQ file.
```
**Example**  
```
    - To analyze a FASTQ file and generate frequency table and plot:  
    
        python script.py -i example.fastq -n 85 -t output_table.csv -p plot.png  
        
    - Specify bases at specific positions simultaneously, for example, from 8 to 16:  
    
        python script.py -i example.fastq -n 85 -t output_table.csv -p plot.png -s 8 16 -sp 8_16_base.png  
        
    - To generate plots directly from an existing frequency table:  
    
        python script.py -a output_table.csv -p plot.png -s 10 15 -sp custom_name.png
```

### 2.  `star_format.py`
Description:  
Parses the .final.out file and extracts the relevant information from the STAR alignment result.

**Usage:**
```bash
python star_format.py DirectoryPath OutputName
```
