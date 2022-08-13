#!/bin/bash

# The script takes three arguments: the filename_1 filename_2 for which the computel 
# should run and the full path to computel/src directory

# The script creates a folder for the computel output in the
# same directory where the files are. That is if the files are in 
# /path/to/file directory and are called ERR242514_1.fastq.gz
#  ERR242514_2.fastq.gz, then the output folder is /path/to/file/ERR242514

#validate.options.R: change somth with r1
# if r1 is single then ....
# if paired then vector
#2 seed parameters -- one single, abother one paired-end

#to remove picard from computel_wrapper.sh
#pipeline.R: How to warn the user of installing packages (not neccessary), create a separate script for paired-end mode
#functions.R: what should I change: base coverage mean for 2. telomere coverages divide by the read_length and then sum 


filename_1=$1
filename_2=$2
file_base=`basename $1`
file_read=`echo "${file_base%%.*}"`

echo file_read is $file_read

computel_dir=$3

if [[ "$computel_dir" == */ ]]
then
    echo $computel_dir has backslash in the end
    computel_dir=${computel_dir%$"/"}
fi


# read_length calculation (Daria's version)
# read_length function, and then call it for 2 filenames
rl_calculation () { 
    lines_count=`zcat $1 | wc -l`
    reads_count=$((lines_count/4))

    if (( reads_count > 10000 )); then
        read_length_total=`zcat $1 | head -40000 | awk 'NR%4==2' | wc -c`
        ((read_length_total -= 10000))
        read_length="$((read_length_total/10000))"
    else
        read_length_total=`zcat $1 | awk 'NR%4==2' | wc -c`
        ((read_length_total -= reads_count))
        read_length="$((read_length_total/reads_count))"
    fi
}
rl_calculation $filename_1
read_length_1=$read_length
echo "read_length_1 is $read_length_1"
rl_calculation $filename_2
read_length_2=$read_length
echo "read_length_2 is $read_length_2"

#read_length function, and then call it for 2 filenames
#output_dir=`echo ${filename%$".fastq"} | rev | cut -d '_' -f2- | rev`
output_dir=$4

echo output_dir is $output_dir


cat > config_file  
scripts.dir	$computel_dir/src/scripts
bowtie.build.path	$computel_dir/bowtie2-2.1.0-linux/bowtie2-build
bowtie.align.path	$computel_dir/bowtie2-2.1.0-linux/bowtie2-align
samtools.path	$computel_dir/samtools-0.1.19-linux/samtools

#read_lengths for each filename
single	T
files.with.prefix	F
file.compression	gz
fastq	$1,$2
read.length $read_length_1,$read_length_2


pattern	TTAGGG
num.haploid.chr	23
min.seed	12
mode.local	F

compute.base.cov	F
##if files are compressed and the base coverage needs to be estimated during unzipping, set "estimate.base.cov" to T
estimate.base.cov	T
genome.length	3244610000

output.dir	$output_dir	

num.proc	4
ignore.err	F
EOF

# echo file_base is $file_base
# echo file_read is $file_read

# echo read_number is $read_number
# echo read_length is $read_length
# echo base_cov is $base_cov
# echo outputdir/file_read is $output_dir/$file_read

if [ ! -d "$output_dir" ]
then
    mkdir $output_dir
    echo $output_dir is created
else
    echo $output_dir exists
fi


# if [ ! -d "$output_dir/$file_read" ]
# then
#    mkdir $output_dir/$file_read
#    echo $output_dir/$file_read is created
# else
#    echo $output_dir/$file_read exists
# fi

Rscript $computel_dir/src/scripts/computel.cmd.R config_file > $output_dir/telseq.log

cp config_file $output_dir                                                                                                                              
