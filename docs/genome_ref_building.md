# BUILDING NEW GENOME REFERENCE

  ## Creating a directory for genome files storage
  ```
  mkdir genome
  mkdir genome/fa
  cd genome/fa
  ```
  
  ## Downloading a genome
  I wasn't able to download a genome directly from the website:
  ```
  wget https://api.ncbi.nlm.nih.gov/datasets/v1/genome/download?filename=GCF_009914755.1.zip&ncbi_phid=322C5C5BD8879F5500003B327174E448.1.m_3.062.1
  ```
  
  The output:
  ```
  --2022-06-08 19:06:16--  https://api.ncbi.nlm.nih.gov/datasets/v1/genome/download?filename=GCF_009914755.1.zip
  Resolving api.ncbi.nlm.nih.gov (api.ncbi.nlm.nih.gov)... 130.14.29.110, 2607:f220:41e:4290::110
  Connecting to api.ncbi.nlm.nih.gov (api.ncbi.nlm.nih.gov)|130.14.29.110|:443... connected.
  HTTP request sent, awaiting response... 405 Method Not Allowed
  2022-06-08 19:06:18 ERROR 405: Method Not Allowed.
  ```
  
  So, I first downloaded FASTA file on my computer and then used an scp command to upload it on server:
  ```
  scp -i ~/.ssh/your_key.rsa ~/Desktop/GCF_009914755.1.zip daria@185.127.66.93:/storage2/users/daria/computel/genome/fa/
  ```

  Then unzipped the file:
  ```
  unzip GCF_009914755.1.zip
  ```
  
  The `/genome/fa/ncbi_dataset/data/GCF_009914755.1` directory was created. It has a FASTA file `GCF_009914755.1_T2T-CHM13v2.0_genomic.fna` that we will further use for genome reference building.
  
  ## Building genome index with bowtie2-build
 
  Creating a directory, where we will store our index:
  ```mkdir genome/index```
  
  Store the index in the genome/index directory with name "t2t-chm13v2.0": 
  ```
  genome/index/t2t-chm13v2.0
  ```
  Command to use bowtie2-build:
  ```
  bowtie2-build FASTAFILEPATH genome/index/t2t-chm13v2.0
  ```
  Before running scripts, create a "log" folder where will be the output information:
  ```mkdir GCF_009914755.1/log```
  
  Building genome reference is a resource-consuming process, so it's important to run it with a slurm resource manager.
  
  Usage steps:
  1) Creation of a .sh file:
      ```touch GenomeReferenceBuilding_extended.sh```
  2) Put the script in this file:
  
      To open this file:
      ```vi GenomeReferenceBuilding_extended.sh```
      
      To start modifying it type `i` (insert mode)
      
      Then you can start typing in this file.
      
      To exit insert mode press esc.
      
      To save the file type `:x`
      
      If you didn't make any changes type `:q` to exit vim:
      
     In the file should be a script that will help us run our resource-consuming programs through the slurm resource manager. Below the sript for slurm we should put our commands to run.
     This is my GenomeReferenceBuilding_extended.sh file that I used for genome reference building with a slurm resource manager:
     
      ```
      #!/bin/bash
      #SBATCH --job-name="ReferenceBuilding"
      #SBATCH --time=10:00:00
      #SBATCH --mem=100GB
      #SBATCH --mail-user=dlaricheva1@gmail.com
      #SBATCH --output=log/ReferenceBuilding2.log
      #SBATCH --mail-type=ALL
      
      srun bowtie2-build /storage2/users/daria/computel/genome/fa/ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna /storage2/users/daria/computel/genome/index/t2t-chm13v2.0
      ```
      NOTE: The script will not run if you are in the wrong folder. You should be in a directory, where your .sh file is, as well as `log` file.
  3)  To run the script:
      ```
      sbatch GenomeReferenceBuilding_extended.sh
      ```
  4)  Command `squeue` helps us see which jobs are currrently running.

  To see the output go to `log` folder and you will see ReferenceBuilding2.log file with the output.

  # Computel usage with the option to align on the new reference
  General usage: 
  ```
  ./computel.sh [options] -1 <fq1> -2 <fq2> -o <outputpath>
  ```
  Option to align on the new reference:
  <-ref> Reference genome index prefix path (generated with bowtie2-build) that you would like to use for more accurate estimation of base coverage.
  
  Command to run computel with a new reference:
  ```
  ./computel.sh -ref /storage2/users/daria/computel/genome/index/t2t-chm13v2.0 -1 src/examples/tel_reads1.fq.gz -2 src/examples/tel_reads2.fq.gz -o mytest_new_ref
  ```
  I performed the same test run but with a new reference.
  
  For comparison of tel.align.bam.coverage.txt files (old test run and test run with a new reference) go to `/storage2/users/daria/computel/mytest_new_ref/align` and to `/storage2/users/daria/computel/mytest/align`