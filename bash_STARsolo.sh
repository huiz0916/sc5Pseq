#!/bin/bash

#base options
while getopts f:r:o:g:a:s: flag
do
        case "$flag" in
                f) forward=${OPTARG};;
                r) reverse=${OPTARG};;
                o) output=${OPTARG};;
		g) genome=${OPTARG};;
		a) annotation=${OPTARG};;
                s) solo_out_prefix=${OPTARG};;
        esac
done

#PROJECT="/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/hui/sc5Pseq"
PROJECT="/cfs/klemming/scratch/h/huizhou/SP_S96_SDHis"
scratch_path="/cfs/klemming/scratch/h/huizhou/SP_S96_SDHis"
#genome_file="/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/hui/genome/S_cerev_S288C/GCF_000146045.2_R64_genomic.fna"
#gtf_file="/cfs/klemming/projects/supr/sllstore2017018/nobackup/projects/hui/genome/S_cerev_S288C/GCF_000146045.2_R64_genomic.gtf"
script_path="/cfs/klemming/home/h/huizhou/script/sc5Pseq"

echo ${forward}
echo ${reverse}
echo ${output}
echo ${genome}
echo ${annotation}
echo ${solo_out_prefix}
#STARsolo

#v1 settings --soloCBposition --soloCBposition 0_8_0_15 0_32_0_39 0_56_0_63 without adapter
#v2 settings --soloCBposition with the adapterseq

#v2
STAR --runThreadN 8 \
        --genomeDir ${genome} \
        --readFilesIn ${forward} ${reverse} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --runDirPerm All_RWX \
         --outFileNamePrefix ${output} \
         --soloType CB_UMI_Complex \
         --soloAdapterSequence AGGTCCTTGGCTTCGC \
         --soloCBposition 2_-9_2_-1 2_16_2_24 2_41_2_49 \
         --soloUMIposition 2_-17_2_-10 \
         --soloAdapterMismatchesNmax 2 \
         --soloCBwhitelist $script_path/metafiles/BC1_v2.txt $script_path/metafiles/BC2_v2.txt $script_path/metafiles/BC3_v2.txt \
         --soloCBmatchWLtype 1MM \
         --soloUMIdedup 1MM_Directional_UMItools \
         --soloBarcodeReadLength 0 \
         --soloFeatures Gene GeneFull \
         --soloCellReadStats Standard \
         --soloOutFileNames ${solo_out_prefix} features.tsv barcodes.tsv matrix.mtx \
         --limitBAMsortRAM 182536110080 \
         --soloMultiMappers Uniform \
         --alignIntronMax 1000 
         #--sjdbOverhang 24
         #--outTmpDir $PROJECT/TMP/STAR_tmp \
         #--sjdbGTFfile ${annotation} \

##paper
:<<!
STAR --runThreadN 24 \
        --genomeDir ${genome} \
        --readFilesIn ${forward} ${reverse}  \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --sjdbGTFfile ${annotation} \
        --outFileNamePrefix ${output} \
        --soloType CB_UMI_Complex \
        --soloAdapterSequence  GTGGCCGATGTTTCGCATCGGCGTACGACT \
        --soloCBposition 2_-8_2_-1 3_1_3_8 3_31_3_38 \
        --soloUMIposition 2_-18_2_-9 \
        --soloCBwhitelist $script_path/metafiles/Round3_barcodes_new4.txt $script_path/metafiles/Round2_barcodes_new4.txt $script_path/metafiles/Round1_barcodes_new5.txt \
        --soloCBmatchWLtype 1MM \
        --soloMultiMappers Uniform \
        --soloFeatures Gene GeneFull \
        --runDirPerm All_RWX \
        --limitBAMsortRAM 9930854492000

#28460880M
#48612416M
:<<!
#v1_3end
STAR --runThreadN 12 \
         --genomeDir ${genome} \
         --readFilesIn ${forward} ${reverse} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --sjdbGTFfile ${annotation} \
         --runDirPerm All_RWX \
         --outFileNamePrefix ${output} 16\
         --soloType CB_UMI_Complex \
         --soloCBposition 0_8_0_15 0_32_0_39 0_56_0_63 \
         --soloUMIposition 0_0_0_7 \
         --soloCBwhitelist $script_path/metafiles/BC1_v1.txt $script_path/metafiles/BC2_v1.txt $script_path/metafiles/BC3_v1.txt \
         --soloCBmatchWLtype EditDist_2 \
         --soloBarcodeMate 0 \
         --soloBarcodeReadLength 0 \
         --soloFeatures Gene GeneFull SJ Velocyto\
         --soloCellReadStats Standard \
         --soloOutFileNames ./solo features.tsv barcodes.tsv matrix.mtx \
         --outTmpDir $PROJECT/STAR_tmp \
         --limitBAMsortRAM 230854492160 \
         --soloMultiMappers Unique

#182536110080
#230854492160

#v1_5end
STAR --runThreadN 12 \
        --genomeDir ${genome} \
        --readFilesIn ${forward} ${reverse} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --sjdbGTFfile ${annotation} \
         --runDirPerm All_RWX \
         --outFileNamePrefix ${output} \
         --soloType CB_UMI_Complex \
         --soloCBposition 0_8_0_15 0_32_0_39 0_56_0_63 \
         --soloUMIposition 0_0_0_7 \
         --soloCBwhitelist $script_path/metafiles/BC1_v1.txt $script_path/metafiles/BC2_v1.txt $script_path/metafiles/BC3_v1.txt \
         --soloCBmatchWLtype EditDist_2 \
         --soloBarcodeMate 0 \
         --soloBarcodeReadLength 0 \
         --soloFeatures Gene GeneFull SJ Velocyto\
         --soloCellReadStats Standard \
         --soloOutFileNames ${solo_out_prefix} features.tsv barcodes.tsv matrix.mtx \
         --clipAdapterType CellRanger4 \
         --outTmpDir $PROJECT/STAR_tmp \
         --limitBAMsortRAM 230854492160 \
         --soloMultiMappers Unique

#2
STAR --runThreadN 12 \
        --genomeDir ${genome} \
        --readFilesIn ${forward} ${reverse} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --sjdbGTFfile ${annotation} \
         --runDirPerm All_RWX \
         --outFileNamePrefix ${output} \
         --soloType CB_UMI_Complex \
         --soloAdapterSequence AGGTCCTTGGCTTCGC \
         --soloCBposition 2_-8_2_-1 2_16_2_23 2_40_2_47 \
         --soloUMIposition 2_-16_2_-9 \
         --soloAdapterMismatchesNmax 3 \
         --soloCBwhitelist $script_path/metafiles/BC1_v1.txt $script_path/metafiles/BC2_v1.txt $script_path/metafiles/BC3_v1.txt \
         --soloCBmatchWLtype 1MM \
         --soloBarcodeMate 0 \
         --soloBarcodeReadLength 0 \
         --soloFeatures Gene GeneFull SJ Velocyto\
         --soloCellReadStats Standard \
         --soloOutFileNames ${solo_out_prefix} features.tsv barcodes.tsv matrix.mtx \
         --outTmpDir $PROJECT/TMP/STAR_tmp \
         --limitBAMsortRAM 182536110080 \
         --soloMultiMappers Unique

STAR --runThreadN 12 \
        --genomeDir ${genome} \
        --readFilesIn ${forward} ${reverse} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --sjdbGTFfile ${annotation} \
         --runDirPerm All_RWX \
         --outFileNamePrefix ${output} \
         --soloType CB_UMI_Complex \
         --soloCBposition 0_56_0_63 \
         --soloUMIposition 0_0_0_7 \
         --soloCBwhitelist $script_path/metafiles/BC3_v1.txt \
         --soloCBmatchWLtype EditDist_2 \
         --soloBarcodeMate 0 \
         --soloBarcodeReadLength 0 \
         --soloFeatures Gene GeneFull SJ Velocyto\
         --soloCellReadStats Standard \
         --soloOutFileNames ${solo_out_prefix} features.tsv barcodes.tsv matrix.mtx \
         --outTmpDir $PROJECT/STAR_tmp \
         --limitBAMsortRAM 230854492160 \
         --soloMultiMappers Unique

STAR --runThreadN 12 \
        --genomeDir ${genome} \
        --readFilesIn ${forward} ${reverse} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
         --sjdbGTFfile ${annotation} \
         --runDirPerm All_RWX \
         --outFileNamePrefix ${output} \
         --soloType CB_UMI_Complex \
         --soloCBposition 0_8_0_15 0_32_0_39 0_56_0_63 \
         --soloUMIposition 0_0_0_7 \
         --soloCBmatchWLtype EditDist_2 \
         --soloBarcodeMate 0 \
         --soloBarcodeReadLength 0 \
         --soloFeatures Gene GeneFull SJ Velocyto\
         --soloCellReadStats Standard \
         --soloOutFileNames ${solo_out_prefix} features.tsv barcodes.tsv matrix.mtx \
         --outTmpDir $PROJECT/STAR_tmp \
         --limitBAMsortRAM 230854492160 \
         --soloMultiMappers Unique \
         --soloCBwhitelist None None None
!