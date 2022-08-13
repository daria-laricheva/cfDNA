#!/bin/bash
#SBATCH --job-name="preprocessing_qc_test"
#SBATCH --time=10000:00:00
#SBATCH --mem=100GB
#SBATCH --mail-user=dlaricheva1@gmail.com
#SBATCH --output=log/preprocessing_qc_nareks_data.log
#SBATCH --mail-type=ALL

# Define an array of samples that have issues with quality:
## SRR119...
fastq_dir="/storage/users/narek/capstone/fastq"
out_dir_119="/storage2/users/daria/cutadapt_119"

bad_quality_samples_119=("SRR1197554_2.fastq.gz" "SRR1197555_2.fastq.gz" "SRR1197556_2.fastq.gz" "SRR1197557_2.fastq.gz" "SRR1197558_2.fastq.gz")
for file1 in ${fastq_dir[@]}/SRR119*_2.fastq.gz; do
    file_name=$(basename ${file1})
    if [[ "${bad_quality_samples_119[@]}" =~ "${file_name}" ]]; then
        echo "file1 = $file1"
        file2="${fastq_dir}/${file_name%%_*}_1.fastq.gz"
        echo "file2 = $file2"
        out1="${out_dir_119}/trimmed_$(basename ${file1})"
        echo "out1 = $out1"
        out2="${out_dir_119}/trimmed_$(basename ${file2})"
        echo "out2 = $out2"
        cutadapt -q 25,25 -m 20 -o ${out1} -p ${out2} ${file1} ${file2}
    fi
done
# fastqc command to check the quality of the output data after trimming.
fastqc ${out_dir_119}/*fastq.gz -o /storage2/users/daria/fastqc_119/
# multiqc command to merge fastqc reports into one.
multiqc /storage2/users/daria/fastqc_119/*zip -o /storage2/users/daria/multiqc_119/

## SRR744...
bad_quality_samples_744=("SRR7442635_2.fastq.gz" "SRR7442636_2.fastq.gz" "SRR7442637_2.fastq.gz" "SRR7442638_2.fastq.gz" "SRR7442639_2.fastq.gz" "SRR7442640_2.fastq.gz" "SRR7442641_2.fastq.gz" "SRR7442642_2.fastq.gz" "SRR7442643_2.fastq.gz" "SRR7442644_2.fastq.gz" "SRR7442645_2.fastq.gz" "SRR7442646_2.fastq.gz" "SRR7442647_2.fastq.gz" "SRR7442648_2.fastq.gz" "SRR7442649_2.fastq.gz" "SRR7442650_2.fastq.gz" "SRR7442651_2.fastq.gz" "SRR7442652_2.fastq.gz" "SRR7442653_2.fastq.gz" "SRR7442654_2.fastq.gz")    
out_dir_744="/storage2/users/daria/cutadapt_744"
# This is the script for paired-end data preprocessing.
for file1 in ${fastq_dir[@]}/SRR744*_2.fastq.gz; do
    # Get the name of SRR*_2.fastq.gz files
    file1_name=$(basename ${file1})
    # Check if the samples from bad_quality_samples array are in a directory with files.
    if [[ "${bad_quality_samples_744[@]}" =~ "${file1_name}" ]]; then 
        echo "file1 = $file1"
        echo
        # Define _1 files based on _2 files' names (file_name%%_* means that we are keeping everything before the underscore (i.e. file's name)).
        file2="${fastq_dir}/${file1_name%%_*}_1.fastq.gz"
        echo "file2 = $file2"
        echo
        out1="${out_dir_744}/trimmed_$(basename ${file1})"
        echo "out1 = $out1"
        echo
        out2="${out_dir_744}/trimmed_$(basename ${file2})"
        echo "out2 = $out2"
        echo
        cutadapt -u 8 -q 20,20 -a AGATCGGAAGAG -A AGATCGGAAGAG -a CTGTCTCTTATA -A CTGTCTCTTATA -n 4 -m 20 -o ${out1} -p ${out2} ${file1} ${file2}
    fi
done 
# fastqc command to check the quality of the output data after trimming.
fastqc ${out_dir_744}/*fastq.gz -o /storage2/users/daria/fastqc_744/
# multiqc command to merge fastqc reports into one.
multiqc /storage2/users/daria/fastqc_744/*zip -o /storage2/users/daria/multiqc_744/

## SRR117...
bad_quality_samples_117=("SRR11742899_1.fastq.gz" "SRR11742899_2.fastq.gz" "SRR11742858_1.fastq.gz" "SRR11742858_2.fastq.gz" "SRR11742869_1.fastq.gz" "SRR11742869_2.fastq.gz" "SRR11742884_1.fastq.gz" "SRR11742884_2.fastq.gz" "SRR11742893_1.fastq.gz" "SRR11742893_2.fastq.gz" "SRR11742911_1.fastq.gz")
out_dir_117="/storage2/users/daria/cutadapt_117"
## Create an array of samples that are already processed.
processed_samples=()
for file_name in ${bad_quality_samples_117[@]}; do
    # Skip files that were already processed.
    if [[ "${processed_samples[@]}" =~ "${file_name%%_*}" ]]; then
        echo "skipping $file_name"
        continue
    fi
    # Define a file path for _1.fastq.gz files. file_name%%_* means that we are keeping everything before the underscore (i.e. file's name).
    file1="${fastq_dir}/${file_name%%_*}_1.fastq.gz"
    echo "file1 = $file1"
    file2="${fastq_dir}/${file_name%%_*}_2.fastq.gz"
    echo "file2 = $file2"
    # If a pair exists then we handle it as paired-end sample.
    if [ -f "$file1" ] && [ -f "$file2" ]; then
        echo "pair exists"
        out1="${out_dir_117}/trimmed_$(basename ${file1})"
        echo "out1 = $out1"
        out2="${out_dir_117}/trimmed_$(basename ${file2})"
        echo "out2 = $out2"
        cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -m 20 -o ${out1} -p ${out2} ${file1} ${file2}
    else
        # If a pair doesn't exist then we handle the sample as a single-end.
        echo "pair doesn't exist"
        # cutadapt for single-end
    fi
    # Keep track of files that were processed to avoid double processing.
    processed_samples+=("${file_name%%_*}")
    echo "processed samples = ${processed_samples[@]}"
done
# fastqc command to check the quality of the output data after trimming.
fastqc ${out_dir_117}/*fastq.gz -o /storage2/users/daria/fastqc_117/
# multiqc command to merge fastqc reports into one.
multiqc /storage2/users/daria/fastqc_117/*zip -o /storage2/users/daria/multiqc_117/