#!/bin/bash
#SBATCH --job-name="computel_117"
#SBATCH --time=1000:00:00
#SBATCH --mem=100GB
#SBATCH --mail-user=dlaricheva1@gmail.com
#SBATCH --output=log/computel_117.log
#SBATCH --mail-type=ALL

trimmed_117="/storage2/users/daria/cutadapt_out_117/*.fastq.gz"
untrimmed_117="/storage/users/narek/capstone/fastq/SRR117*.fastq.gz"
computel_dir="/storage2/users/daria/computel"
users_dir="/storage2/users/daria"
fastq_dir="/storage/users/narek/capstone/fastq"
processed_samples=()
cd /storage2/users/daria/computel

for file in ${trimmed_117[@]}; do
    file_name=$(basename ${file}) 
    echo "files = $file_name"
    file_name_0="${file_name##*trimmed_}"
    if [[ "${processed_samples[@]}" =~ "${file_name_0%%_*}" ]]; then
        echo "skipping $file_name"
        continue
    fi
    file_name_1="${users_dir}/cutadapt_out_117/trimmed_${file_name_0%%_*}_1.fastq.gz"
    file_name_2="${users_dir}/cutadapt_out_117/trimmed_${file_name_0%%_*}_2.fastq.gz"
    echo "file_name_1 = $file_name_1"
    echo "file_name_2 = $file_name_2"
    out="${computel_dir}/computel_out_117_trimmed/${file_name_0%%_*}"
    echo "out = $out"
    ./computel.sh -ref ${computel_dir}/genome/index/t2t-chm13v2.0 -1 ${file_name_1} -2 ${file_name_2} -o ${out}
    processed_samples+=("${file_name_0%%_*}")
    echo "processed samples = ${processed_samples[@]}"
done

processed_samples_untrimmed=()
for untrimmed_file in ${untrimmed_117[@]}; do
    untrimmed_file_name=$(basename ${untrimmed_file})
    echo "untrimmed_files = $untrimmed_file_name"
    file_name_0="${untrimmed_file_name%%_*}"
    if [[ "${processed_samples[@]}" =~ "${file_name_0}" ]]; then
        echo "skipping $untrimmed_file_name"
        continue
    fi
    
    if [[ "${processed_samples_untrimmed[@]}" =~ "${file_name_0}" ]]; then
        echo "skipping $untrimmed_file_name"
        continue
    fi
    file_name_1="${fastq_dir}/${file_name_0}_1.fastq.gz"
    file_name_2="${fastq_dir}/${file_name_0}_2.fastq.gz"
    echo "file_name_1 = $file_name_1"
    echo "file_name_2 = $file_name_2"

    if [ -f "$file_name_1" ] && [ -f "$file_name_2" ]; then
        echo "pair exists"
        out_paired="${computel_dir}/computel_out_117_untrimmed/${file_name_0%%_*}"
        ./computel.sh -ref ${computel_dir}/genome/index/t2t-chm13v2.0 -1 ${file_name_1} -2 ${file_name_2} -o ${out_paired}
    else
        # If a pair doesn't exist then we handle the sample as a single-end.
        echo "pair doesn't exist"
        out_single="${computel_dir}/computel_out_117_single_end_untrimmed/${file_name_0%%_*}"
        ./computel.sh -ref ${computel_dir}/genome/index/t2t-chm13v2.0 -1 ${file_name_1} -o ${out_single}
    fi
    # Keep track of files that were processed to avoid double processing.
    processed_samples_untrimmed+=("${file_name_0}")
    echo "processed_samples_untrimmed = ${processed_samples_untrimmed[@]}"
done
