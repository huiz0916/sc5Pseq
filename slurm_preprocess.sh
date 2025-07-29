#!/bin/bash -l

#SBATCH -A naiss2025-22-734 #
#SBATCH -J qc_starsolo
#SBATCH -p shared #memory shared
#SBATCH -n 1 # Number of tasks
#SBATCH -c 4 # Number of cpus per task
#SBATCH -t 2:00:00

# output log
#SBATCH -o qc_SP_cutadpt.out
#SBATCH -e qc_SP_cutadpt.err

#SBATCH --mail-type ALL --mail-user hui.zhou@scilifelab.se

###export OMP_NUM_THREADS=12
###SBATCH --mem=120G

proj_dir="/cfs/klemming/scratch/h/huizhou/SP2_S96-SK1_SDHis_T3"
fastq_dir=$proj_dir/01_rawdata
trim_dir=$proj_dir/02_cleandata
genome_fsa="/crex/proj/sllstore2017018/nobackup/projects/hui/genome/S_cerev_embl_74/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.fa"
genome_gff="/crex/proj/sllstore2017018/nobackup/projects/hui/genome/S_cerev_embl_74/Saccharomyces_cerevisiae.EF4.74.gtf"

module load bioinfo-tools;wait
module load FastQC;wait
module load MultiQC;wait
module load cutadapt;wait

########### fastqc_raw ###########
echo "########### fastqc raw ###########"
fastqc_dir=$trim_dir/fastqc_raw
if [ ! -d $fastqc_dir ]; then 
    mkdir -p $fastqc_dir
else
    echo ""
    echo "WARNING: the contents of the directory $fastqc_dir may be overwritten." 
    echo ""
fi

fastqcdir=`which fastqc`
if [ -z $fastqcdir ]; then
	echo "Could not find program FastQC. Make sure you have it installed or loaded."
	exit -1
else
	echo ""
fi

for f in $fastq_dir/*.gz
do
	echo "$f: fastqc"
	fastqc -t 4 -o $fastqc_dir $f
	wait
done
echo "FastQC finished for all fastq files in $fastq_dir. Output written to $trim_dir"


########### multiqc raw ###########
echo "########### multiqc raw ###########"

multiqc_dir=$trim_dir/multiqc_raw

if [ ! -d $multiqc_dir ]; then 
mkdir $multiqc_dir
else
    echo ""
    echo "WARNING: the contents of the directory $multiqc_dir may be overwritten." 
    echo ""
fi

multiqcdir=`which multiqc`
if [ -z $multiqcdir ]; then
    echo "Could not find program MULTIQC. Make sure you have it installed or loaded."
    exit -1
else
    echo ""
fi

multiqc -d -o $multiqc_dir $fastqc_dir

echo "MULTIQC finished for all fastqc files in $fastqc_dir"

########### cutadapt ###########
echo "########### cutadapt ###########"
moduledir=`which cutadapt`
if [ -z $moduledir ]; then
    echo "Could not find program cutadapt. Make sure you have it installed or loaded."
    exit -1
else
    echo ""
fi

cutadapt_dir=$trim_dir/fastq_trimmed
if [ ! -d $cutadapt_dir ]; then 
    mkdir $cutadapt_dir
else
    echo "WARNING: the contents of the directory $cutadapt_dir may be overwritten."
fi

echo ""

for f in $fastq_dir/*_R1_001.fastq.gz
do
    fpath=${f%/*}
    r1=${f##*/}
    r2=${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}

    if [ -f "$fpath/$r2" ]; then
        trimed_r1=${r1/.fastq.gz/.trim.fastq.gz}
        trimed_r2=${r2/.fastq.gz/.trim.fastq.gz}
        echo ""
        echo "$r1,$r2: trimming"
        #for the low quality read 2 of first 13 bases
        #cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -U 13 -Q 20,0 --minimum-length 73:25 -e 0.2 --nextseq-trim=20 -o $cutadapt_dir/$trimed_r1 -p $cutadapt_dir/$trimed_r2 $fpath/$r1 $fpath/$r2
        cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT --minimum-length 73:35 -e 0.2 --nextseq-trim=20 -o $cutadapt_dir/$trimed_r1 -p $cutadapt_dir/$trimed_r2 $fpath/$r1 $fpath/$r2 -j 4

        wait
        success=$?

        if [ $success -eq 0 ]; then
            echo "$r1,$r2: successfully trimmed"
            echo ""
        else
            echo "Cutadapt finished with non-zero exit status $success for file $f"
            exit $success
        fi
    else
        echo "Matching R2 file not found for $fpath/$r1"
    fi
done

echo "Successfully finished adapter trimming for all files"
echo ""
#fastq_dir=$cutadapt_dir

:<<!
### this is for 5 end with TSO sequence
echo "########### cutadapt 5' TSO ###########"
for f in $fastq_dir/*_*5E*R1_001.trim.fastq.gz
do
    fpath=${f%/*}
    r1=${f##*/}
    r2=${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}

    if [ -f "$fpath/$r2" ]; then
        trimed_r1=${r1/.fastq.gz/.trim_TSO.fastq.gz}
        trimed_r2=${r2/.fastq.gz/.trim_TSO.fastq.gz}
        echo ""
        echo "$r1,$r2: trimming"
        cutadapt -G "^WWGG" --minimum-length 73:30 -e 0.2 --nextseq-trim=20 -o $cutadapt_dir/$trimed_r1 -p $cutadapt_dir/$trimed_r2 $fpath/$r1 $fpath/$r2 -j 4
        
        wait
        success=$?

        if [ $success -eq 0 ]; then
            echo "$r1,$r2: successfully trimmed"
            echo ""
        else
            echo "Cutadapt finished with non-zero exit status $success for file $f"
            exit $success
        fi
    else
        echo "Matching R2 file not found for $fpath/$r1"
    fi
done

echo "Successfully finished adapter trimming for all files"
echo ""
fastq_dir=$cutadapt_dir
fastq_dir=$trim_dir/clean_data
!

fastq_dir=$trim_dir/fastq_trimmed
########### fastqc ###########
echo "########### fastqc ###########"

module load bioinfo-tools FastQC
wait

fastqc_dir=$trim_dir/fastqc_trim_3E

echo "INFO: fastqc input" 
echo "  fastq: $fastq_dir"
echo "  output: $fastqc_dir"
echo ""

if [ ! -d $fastqc_dir ]; then 
    mkdir $fastqc_dir
else
    echo ""
    echo "WARNING: the contents of the directory $fastqc_dir may be overwritten."
    echo ""
fi


for f in $fastq_dir/*.gz
do
    echo "$f: fastqc"
    fastqc -t 4 -o $fastqc_dir $f
    wait
done
echo "FastQC finished for all fastq files in $fastq_dir. Output written to $trim_dir"

########### multiqc ###########
echo "########### multiqc ###########"

module load bioinfo-tools MultiQC
wait

multiqc_dir=$trim_dir/multiqc_trim_3E

echo "INFO: fastqc input" 
echo "  fastqc: $fastqc_dir"
echo "  output: $multiqc_dir"
echo ""


if [ ! -d $multiqc_dir ]; then 
    mkdir $multiqc_dir
else
    echo ""
    echo "WARNING: the contents of the directory $multiqc_dir may be overwritten."
    echo ""
fi

multiqcdir=`which multiqc`
if [ -z $multiqcdir ]; then
    echo "Could not find program MULTIQC. Make sure you have it installed or loaded."
else
    echo ""
fi

multiqc -d -o $multiqc_dir $fastqc_dir

echo "MULTIQC finished for all fastqc files in $fastqc_dir"
