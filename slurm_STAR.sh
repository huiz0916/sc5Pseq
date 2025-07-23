#!/bin/bash -l
# The -l above is required to get the full environment with modules

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A naiss2025-22-734 #naiss2025-22-440

# The name of the script is myjob
#SBATCH -J staro

# The partition
#SBATCH -p shared #memory shared

# Number of tasks
#SBATCH -n 1

# Number of cpus per task
#SBATCH -c 8

# X hours wall clock time will be given to this job
#SBATCH -t 01:00:00

# output log
#SBATCH -o star.out
#SBATCH -e star.err

#SBATCH --mail-type ALL --mail-user hui.zhou@scilifelab.se

##export OMP_NUM_THREADS=12
module load rnastar/2.7.11b


PROJECT="/cfs/klemming/scratch/h/huizhou/SP2_S96-SK1_SDHis"
fastq_dir="$PROJECT/02_cleandata/fastq_trimmed"
#genome_file="/cfs/klemming/scratch/h/huizhou/SP_lab_diff/08_star_index/"
genome_file="/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/hui/genome/S_cerev_embl_74/"
index_path="/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/hui/genome/S_cerev_embl_74_star_len25_index/"
#gtf_file="/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/hui/genome/S_cerev_S288C/GCF_000146045.2_R64_genomic.gtf"
gtf_file="/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/hui/genome/S_cerev_embl_74/Saccharomyces_cerevisiae.EF4.74.gtf"
script_path="/cfs/klemming/home/h/huizhou/script/sc5Pseq"

#v1

srun $HOME/script/sc5Pseq/bash_STARsolo.sh \
     -f $fastq_dir/S96-SK1_SDHis_1K3E_R1T1_S19_R2_001.trim.fastq.gz \
     -r $fastq_dir/S96-SK1_SDHis_1K3E_R1T1_S19_R1_001.trim.fastq.gz \
     -o $PROJECT/SP2_S96-SK1_3E_dedup \
     -s SP2_3E \
     -g $genome_file \
     -a $gtf_file


:<<!
for f in $PROJECT/r2/04_trimmed_reads_r2_tso_only/B*R1*.fastq.gz;
do
     fpath=${f%/*}
     r1=${f##*/}
     r2=${r1/_R1.trim2_tso_only.fastq.gz/_R2.trim2_tso_only.fastq.gz}
     solo_out_prefix=${r1/_R1.trim2_tso_only.fastq.gz/}
     echo "STARsolo runing: ${solo_out_prefix}"

     srun $HOME/script/sc5Pseq/bash_STARsolo.sh \
     -f $fpath/$r2 \
     -r $fpath/$r1 \
     -o $PROJECT/STARsolo_${solo_out_prefix}_tso_only/${solo_out_prefix}_ \
     -s ${solo_out_prefix}_ \
     -g $genome_file \
     -a $gtf_file
done
!
