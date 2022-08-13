# Quality Control and Data Trimming (Narek's data)
  ## Decription of the data
  In the directory ```/storage/users/narek/capstone/fastq/``` there are three experiments: SRR117...(ncbi link: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA554329), SRR119...(ncbi link: https://www.ncbi.nlm.nih.gov/sra/SRR1197524), SRR744...(ncbi link: https://www.ncbi.nlm.nih.gov/sra/SRR7442635).
  ## Initial quality control of the datasets
  The command for quality check: ``` fastqc /storage/users/narek/capstone/fastq/*fastq.gz -o /storage2/users/daria/fastqc_narek_data_out/ ```.
  
  MultiQC can be used in order to merge all the reports into one. The command to merge reports: ``` multiqc /storage2/users/daria/fastqc_narek_data_out/*zip -o /storage2/users/daria/multiqc_narek_data_out/ ```
  
  ### Results
  You may found results of the QC of raw data here: ``` /storage2/users/daria/multiqc_narek_data_out/multiqc_report.html ```
  
  <img width="667" alt="Снимок экрана 2022-07-31 в 12 08 34" src="https://user-images.githubusercontent.com/106603217/182035295-8ee4822e-6857-45d6-b6d8-ac5d0c45e93d.png">
  
  <img width="682" alt="Снимок экрана 2022-07-31 в 13 13 19" src="https://user-images.githubusercontent.com/106603217/182037774-2dfa6257-8bf8-406a-83f1-7e58f4021b96.png">
  
  <img width="666" alt="Снимок экрана 2022-07-31 в 12 09 17" src="https://user-images.githubusercontent.com/106603217/182035325-34b1d29e-ac8e-4d55-8d2e-a60a6aaefe6c.png">

<img width="671" alt="Снимок экрана 2022-07-31 в 13 11 06" src="https://user-images.githubusercontent.com/106603217/182037684-54eb5581-d2a7-4ca1-a6b2-ea42ff0e341e.png">

<img width="675" alt="Снимок экрана 2022-07-31 в 13 11 54" src="https://user-images.githubusercontent.com/106603217/182037712-d1290a4d-0f4e-4c5d-8c32-04fd7c940171.png">

Judging by sequence quality and adapter content, we can see which samples need adapter trimming and quality trimming.

For quality/adapter trimming we will use Cutadapt (https://cutadapt.readthedocs.io/en/stable/guide.html).

To find which adapters are present in samples, I went through individal fastqc reports. The path is: ``` /storage2/users/daria/fastqc_narek_data_out/ ```. "Adapter content" shows us which adapters present in samples (Illumina Universal Adapter, Nextera Transposase Sequence, etc.). Here is the website with adapter sequences that I used to perform adapter trimming: https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/FastQC_Adapter_Kmer_files_fDG.htm.

Concerning SRR744, we can see that some samples have a quality drop, starting from 0 to 8 bp. This means, we can cut out first 8 bp (``` cutadapt -u 8 ```). We will see it further in the script. Also, we consider to keep only reads with read length more than 20 bp (``` cutadapt -m 20 ```) since it affects downstream analysis.

Since there are many samples, we need to write a script to run Cutadapt on all of them.

The script: https://github.com/abi-am/cfDNA/blob/main/src/bash-scripts/cutadapt_nareks_data.sh

Data preprocessing (quality trimming, adapter trimming) takes some time (around 20 hours totally), so we run it with the slurm.

Once it finishes running, we can find our reports here: ``` /storage2/users/daria/multiqc_119/ ```, ``` /storage2/users/daria/multiqc_744/ ```, ``` /storage2/users/daria/multiqc_117/ ```.

Data is preprocessed now.