# M6Allele Pipeline & M6Allele algorithm

## Introduction
We have developed an algorithm called **M6Allele** for identifying allele-specific m6A modifications. To facilitate its usage by researchers, we have also encapsulated our analysis process into a pipeline. You can learn more about the pipeline and the algorithm's usage from the following two modules:
* [Pipeline](#m6allele-pipeline)
* [M6Allele algorithm](#m6allele-10)

## M6Allele Pipeline
### PARAMETER INTRODUCTION
* `-g/--gtf` : required, the file name of your own GTF file
* `-fa/--fasta` : required, the file name of your own reference genome file
* `-sf/--skip_fastqc` : optional, whether to skip the fastqc phase. Default false
* `-se/--single_end` : required, the fastq files are single-end or paired-end sequencing
* `-vg/--varscan_or_gatk` : optional, use VarScan or GATK to call snp, v: VarScan, g: GATK. Default v
* `-f/--function` : required, the function of M6Allele, including `AseGeneDetection`, `AsmPeakDetection`, `SampleSpecificASM`. Please refer to [M6Allele 1.0](#m6allele-10) for specific explanation
* `-s/--sample` : required, the name of the file containing the sample name to be processed
* `-gzip/--is_gzip` : required, whether the fastq file is compressed
* `-db/--dbSnp` : optional, the name of dbSNP vcf file
* `-h/--help` : help message of the pipeline

### USAGE
#### Overview
1. Install [docker(v24.0.7)](https://www.docker.com/get-started/)
2. Download a compressed docker image file from this [link](https://renlab.oss-cn-shenzhen.aliyuncs.com/M6Allele/m6allelepipe.tar.gz) and import it using the following command:
   ```shell
      cd your_compressed_file_directory
      gunzip m6allelepipe.tar.gz
      docker load -i m6allelepipe.tar
   ```
    * If you're unable to down the image from above link, you can download the required files for local packaging images from [here](https://renlab.oss-cn-shenzhen.aliyuncs.com/M6Allele/docker.tar.gz) and then build the image locally
        ```shell
            # build command
            docker build -t your_image_name .
        ```
3. Assuming your current working directory is `your_work_directory`, you need to create the following subdirectories or files and place the required files:
   * `fastq` : **required**, containing the fastq files you need to process
   * `scripts` : **required**, containing three script files: main.sh, main.py and getPeakInMeT.R, which you can download from this [link](https://renlab.oss-cn-shenzhen.aliyuncs.com/M6Allele/scripts.tar.gz)
   * `reference` : **required**, containing references files required for the processing workflow
        * `your_gtf_file.gtf` : we provide our default fasta file, you can download it from [this](https://renlab.oss-cn-shenzhen.aliyuncs.com/M6Allele/gtfAndfasta.tar.gz)
        * `your_reference_genome_file.fa` : we provide our default gtf file, you can download it from [this](https://renlab.oss-cn-shenzhen.aliyuncs.com/M6Allele/gtfAndfasta.tar.gz)
        * If you want to use GATK for SNP calling, you need to provide the dbSNP dataset. Here, we provide our default dbSNP dataset. You can download the `GCF_000001405.39.dbsnp.vcf.gz`from [here](https://renlab.oss-cn-shenzhen.aliyuncs.com/M6Allele/GCF_000001405.39.dbsnp.vcf.gz) and `GCF_000001405.39.dbsnp.vcf.gz.tbi` files from [here](https://renlab.oss-cn-shenzhen.aliyuncs.com/M6Allele/GCF_000001405.39.dbsnp.vcf.gz.tbi) and place them in the `reference` folder. If you want to use a different version of the dbSNP dataset, please follow these steps:
          * Download the GCF_XXX.vcf.gz and GCF_XXX.vcf.gz.tbi files from the [dbSNP database](https://ftp.ncbi.nih.gov/snp/latest_release/VCF/), and download the [chromosome conversion files](https://ftp.ncbi.nih.gov/genomes/all/GCF/000/001/405/) corresponding to the above-mentioned files. Process the downloaded dbSNP files using the following commands: 
              ```shell
                # convert the chromosome names in the GCF_XXX.vcf.gz file to 1, 2, ..., X, Y
                bcftools annotate --rename-chrs chromosome_conversion.txt --threads 10 -Oz -o your_new_process_dbsnp_file.vcf.gz your_downloaded_dbsnp_file.vcf.gz
                # generate .tbi file
                bcftools index -t your_new_process_dbsnp_file.vcf.gz
              ```
          * Move the resulting new dbSNP database `your_new_process_dbsnp_file.vcf.gz` and `your_new_process_dbsnp_file.vcf.gz.tbi` files to the `reference` folder
          * You also need to provide .fai and .dict index files for the reference genome .fa file. If your fasta file was downloaded from the link we provided, you will also download the corresponding .fai and .dict files
          * However, If you have your own fasta file, you can generate the corresponding .fai and .dict files using the following commands. Then move .fa and .dict file to `reference` folder
            ``` shell
                # Please ensure that the .fai, .dict, and .fa files have the same prefix in their filenames
                # .fai generate command
                samtools faidx your_reference_genome_file.fa
                # .dict generate command
                gatk CreateSequenceDictionary -R your_reference_genome_file.fa -O your_reference_genome_file.dict
            ```
   * `your_sample.txt`: **required**, containing the sample names you want to process, which are the prefixes of the fastq files. **Each line is separated by either space or tab.** According to the M6Allele function, there are three formats:
        * `AseGeneDetection` : **Each line represents the name of an RNA-seq sample**. If there are multiple duplicates, use multiple lines to represent them
        * `AsmPeakDetection` : **Each line consists of two columns, representing the INPUT sample name and the corresponding IP sample name for MeRIP-seq**. If there are multiple duplicates, use multiple lines to represent them
        * `SampleSpecificASM` : **Each line consists of four columns, representing the INPUT sample name and the corresponding IP sample name for MeRIP-seq of sample 1, as well as the INPUT sample name and the corresponding IP sample name for MeRIP-seq of sample 2**. If there are multiple duplicates, use multiple lines to represent them

### Specific Example:
Here, we have listed several specific examples of using the pipeline. If you have other requirements, you can achieve them by combining different parameters.
### 1. To use VarScan for calling SNPs and detecting ASE genes

**data dependency:**\
`your_work_directory`: there are the following subfolders and files
* `fastq` : It contains two files:
  * input1.fastq.gz
  * input2.fastq.gz
* `scripts` : containing three script files: main.sh, main.py and getPeakInMeT.R
* `reference` : 
  * your_fasta_file.fa
  * your_gtf_file.gtf
* `sample.txt` : It contains two lines:
  * First line: input1
  * Second line: input2

**example:**
```shell
  docker run -v /path/to/your_work_directory:/data renlab303/m6allelepipe -f AseGeneDetection -s sample.txt -gzip true -se true -fa your_fasta_file.fa -g your_gtf_file.gtf
```

### 2. To use GATK for calling SNPs and detecting ASM m6A signals
**data dependency:**\
`your_work_directory`: there are the following subfolders and files
* `fastq` : It contains eight files: 
  * input1_1.fastq.gz
  * input1_2.fastq.gz
  * ip1_1.fastq.gz
  * ip1_2.fastq.gz 
  * input2_1.fastq.gz
  * input2_2.fastq.gz
  * ip2_1.fastq.gz
  * ip2_2.fastq.gz
* `scripts` : containing three script files: main.sh, main.py and getPeakInMeT.R
* `reference` :
  * your_gtf_file.gtf
  * your_fasta_file.fa
  * your_fasta_file.fa.fai
  * your_fasta_file.dict
  * your_dbSNP_vcf_file.vcf.gz
  * your_dbSNP_vcf_file.vcf.gz.tbi
* `sample.txt` : It contains two lines
  * input1&emsp;ip1
  * input2&emsp;ip2

**example:**
```shell
  docker run -v /path/to/your_work_directory:/data renlab303/m6allelepipe -f AsmPeakDetection -s sample.txt -gzip true -se false -g your_gtf_file.gtf -fa your_fasta_file.fa -vg g -db your_dbSNP_vcf_file.vcf.gz
```

### 3. To use GATK for calling SNPs and detecting Sample-specific ASM m6A signals
**data dependency:**\
`your_work_directory`: there are the following subfolders and files
* `fastq` : It contains eight files: 
  * sample1_input1.fastq.gz
  * sample1_ip1.fastq.gz
  * sample2_input1.fastq.gz
  * sample2_ip1.fastq.gz
  * sample1_input2.fastq.gz
  * sample1_ip2.fastq.gz
  * sample2_input2.fastq.gz
  * sample2_ip2.fastq.gz 
* `scripts` : containing three script files: main.sh, main.py and getPeakInMeT.R
* `reference` : It contains following files:
    * your_gtf_file.gtf
    * your_fasta_file.fa
    * your_fasta_file.fa.fai
    * your_fasta_file.dict
    * your_dbSNP_vcf_file.vcf.gz
    * your_dbSNP_vcf_file.vcf.gz.tbi
* `sample.txt` : It contains two lines
    * sample1_input1&emsp;sample1_ip1&emsp;sample2_input1&emsp;sample2_ip1
    * sample1_input2&emsp;sample1_ip2&emsp;sample2_input2&emsp;sample2_ip2

**example:**
```shell
  docker run -v /path/to/your_work_directory:/data renlab303/m6allelepipe -f SampleSpecificASM -s sample.txt -gzip true -se true -g your_gtf_file.gtf -fa your_fasta_file.fa -vg g -db your_sbSNP_vcf_file.vcf.gz
```

### Pipeline overview

This pipeline is built using shell scripts and integrates tools as follows:

* **Quality control and preprocessing of raw data**
    * [fastp](https://github.com/OpenGene/fastp): quality trimming and adapter clipping
    * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): generate quality reports
* **Build STAR index**
  * [STAR](https://github.com/alexdobin/STAR): build index
* **Read alignment**
    * [STAR](https://github.com/alexdobin/STAR): Spliced Transcripts Alignment to a Reference
    * [Samtools](http://www.htslib.org/): Reads sort and remove duplicates
* **SNP calling**
    * [VarScan](https://varscan.sourceforge.net/): Call SNPs from MeRIP-seq INPUT sample
    * [bcftools](http://www.htslib.org/doc/1.1/bcftools.html): mark the result of VarScan
    * [vcftools](https://vcftools.github.io/): filter the result of bcftools
    * [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels): Call SNPs from MeRIP-seq INPUT sample
* **Peak calling**
    * [MeTPeak](https://github.com/compgenomics/MeTPeak): a novel, graphical model-based peak-calling method
    * [BEDTools](https://bedtools.readthedocs.io/en/latest/): using "mergeBed" function
* **ASE or ASM m6A detection**
    * [M6Allele](#M6Allele 1.0): A toolkit for detection of allele-specific RNA N6-methyladenosine modifications 

## M6Allele 1.0

### HARDWARE/SOFTWARE REQUIREMENTS
* Java 1.8
* Windows / Linux / Mac OS

### INSTALLATION
* clone the repo,
```
git clone https://github.com/Jakob666/allele-specificM6A.git
```
* target JAR package
```
cd ./allele-specificM6A
```
make sure the directory contains `M6Allele.jar`.

### USAGE
#### Overview
#### Tools Introduction
M6Allele.jar provides the following tools:

Tool | Function
---|---
AseGeneDetection|detect allele-specific expression (ASE) genes (one sample test)
AsmPeakDetection|detect allele-specific modification (ASM) m6A signals  (one sample test)
SampleSpecificASM|detect sample-specific ASM m6A signals  (paired sample test)

#### parameters description
1. **AseGeneDetection**
   - `-vcf/--vcf_file` : required, VCF format file generate by RNA-seq or MeRIP-seq data SNP calling process
   - `-g/--gtf` : required, GTF annotation file
   - `-bam/--bam_files` : required, the alignment file of the FASTQ file, if you have multiple BAM files, please separate them with commas
   - `-bai/--bai_files` : required, the index file of the bam file, If you have multiple BAI files, please separate them with commas
   - `-o/--output` : optional, ASE gene test output file, default `./aseGene.txt`
   - `-rc/--reads_coverage` : optional, reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR), default 10
   - `-s/--sampling`: optional, MCMC sampling times, larger than 500, default 50000
   - `-b/--burn` : optional, MCMC burn-in times, more than 100 and less than sampling times. Default 10000
   - `-t/--thread` : optional, thread number for running test. Default 2
   - `-h/--help` : help message of AseGeneDetection

2. **AsmPeakDetection**
   - `-bed/--peak_bed_file` : required, peak calling output result in BED format
   - `-vcf/--vcf_file` : required, VCF format file generate by RNA-seq or MeRIP-seq data SNP calling process
   - `-g/--gtf` : required, GTF annotation file
   - `-inputBam/--input_bam_file` : required, the alignment file of the INPUT sample, if you have multiple BAM files, please separate them with commas
   - `-inputBai/--input_bai_file` : required, the index file of the INPUT bam file, if you have multiple BAI files, please separate them with commas
   - `-ipBam/--ip_bam_file` : required, the alignment file of the Ip sample
   - `-ipBai/--ip_bai_file` : required, the index file of the Ip bam file
   - `-o/--output` : optional, ASM m6A signal test output file, default `./asmPeak.txt`
   - `-rc/--reads_coverage` : optional, reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR), default 10
   - `-s/--sampling`: optional, MCMC sampling times, larger than 500, default 50000
   - `-b/--burn` : optional, MCMC burn-in times, more than 100 and less than sampling times. Default 10000
   - `-t/--thread` : optional, thread number for running test. Default 2
   - `-h/--help` : help message of AsmPeakDetection
   
3. **SampleSpecificASM**
   - `-s1Vcf/--sample1VcfFile` : required, VCF format file generate by sample1 RNA-seq or MeRIP-seq data SNP calling process
   - `-s2Vcf/--sample2VcfFile` : required, VCF format file generate by sample2 RNA-seq or MeRIP-seq data SNP calling process
   - `-bed/--mergePeakBedFile` : required, the merge result of peak calling for Sample 1 and Sample 2 in BED format
   - `-g/--gtf` : required, GTF annotation file
   - `-s1InputBam/--s1InputBamFiles` : required, the alignment file of the sample1 INPUT sample, if you have multiple BAM files, please separate them with commas
   - `-s1InputBai/--s1InputBaiFiles` : required, the index file of the sample1 INPUT sample bam file, if you have multiple BAI files, please separate them with commas
   - `-s2InputBam/--s2InputBamFiles` : required, the alignment file of the sample2 INPUT sample, if you have multiple BAM files, please separate them with commas
   - `-s2InputBai/--s2InputBaiFiles` : required, the index file of the sample2 INPUT sample bam file, if you have multiple BAI files, please separate them with commas
   - `-s1IpBam/--s1IpBamFiles` : required, the alignment file of the sample1 Ip sample, if you have multiple BAM files, please separate them with commas
   - `-s1IpBai/--s1IpBaiFiles` : required, the index file of the sample1 Ip sample bam file, if you have multiple BAI files, please separate them with commas
   - `-s2IpBam/--s2IpBamFiles` : required, the alignment file of the sample2 Ip sample, if you have multiple BAM files, please separate them with commas
   - `-s2IpBai/--s2IpBaiFiles` : required, the index file of the sample2 Ip sample bam file, if you have multiple BAI files, please separate them with commas
   - `-o/--output` : optional, Sample-specific test output directory, default `.`
   - `-rc/--reads_coverage` : optional, reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR), default 10
   - `-s/--sampling`: optional, MCMC sampling times, larger than 500, default 50000
   - `-b/--burn` : optional, MCMC burn-in times, more than 100 and less than sampling times. Default 10000
   - `-t/--thread` : optional, thread number for running test. Default 2
   - `-threshold/--significantThreshold` : optional, the threshold for determining whether it is a sample-specific Asm6A modification, default 0.05
   - `-h/--help` : help message of SampleSpecificASM


### 1. Allele-specific expression (ASE) gene detection (one sample test)
**data dependency**:
1. VCF format file generate by SNV calling process of `RNA-seq data` or `MeRIP-seq INPUT data` (required, the format of the file is described below)
2. GTF file (required)
3. bam file (required)
4. bai file (required, the index file of bam file)

**examples**:\
suppose here exists files below:
1. human genome GTF file `/path/to/Homo_sapiens.GRCh38.93.chr.gtf`
2. VCF format file generate by RNA data `/path/to/rna_filtered.vcf`
3. bam files `/path/to/repeat1.bam (/path/to/repeat2.bam)`
4. bai files `/path/to/repeat1.bam.bai (/path/to/repeat2.bam.bai)`

* detect ASE gene
```
# command
java -jar ./M6Allele.jar AseGeneDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -vcf /path/to/rna_filtered.vcf 
     -bam /path/to/repeat1.bam,/path/to/repeat2.bam
     -bai /path/to/repeat1.bam.bai,/path/to/repeat2.bam.bai
     -o /path/to/output_file 
     -t 6
```


### 2. Allele-specific modification (ASM) m6A signal detection (one sample test)
**data dependency**:
1. GTF format file
2. VCF format file generate by SNV calling process of `RNA-seq data` or `MeRIP-seq INPUT data` (required, the format of the file is described below)
3. BED format peak calling result generate by `MeRIP-seq data` (required, the format of the file is described below)
4. The bam file of `MeRIP-seq INPUT data` (required)
5. The bai file of `MeRIP-seq INPUT data` (required)
6. The bam file of `MeRIP-seq Ip data` (required)
7. The bai file of `MeRIP-seq Ip data` (required)

**examples**:\
suppose here exists files below:
1. human genome GTF file `/path/to/Homo_sapiens.GRCh38.93.chr.gtf`
2. VCF format file generate by RNA data `/path/to/rna_filtered.vcf`
3. BED format file generate by peak calling process `/path/to/peak.bed`
4. Bam files generate by `MeRIP-seq INPUT data` `/path/to/repeat1_input.bam (/path/to/repeat2_input.bam)`
5. Bai files generate by `MeRIP-seq INPUT data` `/path/to/repeat1_input.bam.bai (/path/to/repeat2_input.bam.bai)`
6. Bam files generate by `MeRIP-seq Ip data` `/path/to/repeat1_ip.bam (/path/to/repeat2_ip.bam)`
7. Bai files generate by `MeRIP-seq Ip data` `/path/to/repeat1_ip.bam.bai (/path/to/repeat2_ip.bam.bai)`

* detect ASM m6A signal
```
# command
java -jar ./M6Allele.jar AsmPeakDetection 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -inputBam /path/to/repeat1_input.bam,/path/to/repeat2_input.bam
     -inputBai /path/to/repeat1_input.bam.bai,/path/to/repeat2_input.bam.bai
     -ipBam /path/to/repeat1_ip.bam,/path/to/repeat2_ip.bam
     -ipBai /path/to/repeat1_ip.bam.bai,/path/to/repeat2_ip.bam.bai
     -bed /path/to/peak.bed 
     -vcf /path/to/rna_filtered.vcf 
     -o /path/to/output_file 
     -t 6
```

### 3. Sample-specific ASM m6A signal detection (paired sample test)
**data dependency**:
1. GTF format file
2. paired sample VCF format files generate by SNV calling process of `RNA-seq data` or `MeRIP-seq INPUT data` (required, the format of the file is described below)
3. paired sample BED format peak calling results generate by `MeRIP-seq data` (required, the format of the file is described below)
   * After obtaining the m6A peak calling results for two samples separately, you need to merge the results using bedtools.
4. paired sample Bam files generate by `MeRIP-seq INPUT data`(required)
5. paired sample Bai files generate by `MeRIP-seq INPUT data`(required)
6. paired sample Bam files generate by `MeRIP-seq Ip data`(required)
7. paired sample Bai files generate by `MeRIP-seq Ip data`(required)

**examples**:\
suppose here exists files below:
1. human genome GTF file `/path/to/Homo_sapiens.GRCh38.93.chr.gtf`
2. paired sample VCF format files generate by RNA data `/path/to/sample1_rna_filtered.vcf` & `/path/to/sample2_rna_filtered.vcf`
3. paired sample BED format files generate by peak calling process `/path/to/merge_peak.bed`
4. paired sample Bam files generate by MeRIP-seq INPUT data `/path/to/sample1_repeat1_input.bam (/path/to/sample1_repeat2_input.bam)` & `/path/to/sample2_repeat1_input.bam (/path/to/sample2_repeat2_input.bam)`
5. paired sample Bai files generate by MeRIP-seq INPUT data `/path/to/sample1_repeat1_input.bam.bai (/path/to/sample1_repeat2_input.bam.bai)` & `/path/to/sample2_repeat1_input.bam.bai (/path/to/sample2_repeat2_input.bam.bai)`
6. paired sample Bam files generate by MeRIP-seq Ip data `/path/to/sample1_repeat1_ip.bam (/path/to/sample1_repeat2_ip.bam)` & `/path/to/sample2_repeat1_ip.bam (/path/to/sample2_repeat2_ip.bam)`
7. paired sample Bai files generate by MeRIP-seq Ip data `/path/to/sample1_repeat1_ip.bam.bai (/path/to/sample1_repeat2_ip.bam.bai)` & `/path/to/sample2_repeat1_ip.bam.bai (/path/to/sample2_repeat2_ip.bam.bai)`

* detect sample-specific ASM m6A signal
```
# command
java -jar ./M6Allele.jar SampleSpecificASM 
     -g /path/to/Homo_sapiens.GRCh38.93.chr.gtf 
     -bed /path/to/merge_peak.bed
     -s1Vcf /path/to/sample1_rna_filtered.vcf 
     -s2Vcf /path/to/sample2_rna_filtered.vcf 
     -s1InputBam /path/to/sample1_repeat1_input.bam,/path/to/sample1_repeat2_input.bam
     -s1InputBai /path/to/sample1_repeat1_input.bam.bai,/path/to/sample1_repeat2_input.bam.bai
     -s1IpBam /path/to/sample1_repeat1_ip.bam,/path/to/sample1_repeat2_ip.bam
     -s1IpBai /path/to/sample1_repeat1_ip.bam.bai,/path/to/sample1_repeat2_ip.bam.bai
     -s2InputBam /path/to/sample2_repeat1_input.bam,/path/to/sample2_repeat2_input.bam
     -s2InputBai /path/to/sample2_repeat1_input.bam.bai,/path/to/sample2_repeat2_input.bam.bai
     -s2IpBam /path/to/sample2_repeat1_ip.bam,/path/to/sample2_repeat2_ip.bam
     -s2IpBai /path/to/sample2_repeat1_ip.bam.bai,/path/to/sample2_repeat2_ip.bam.bai
     -o /path/to/output_dir
     -t 6
```

### FORMAT DECLARATION
### 1. VCF generate by SNP calling of RNA and MeRIP sequencing data
At least 2 columns,
* \#CHROM: chromosome number, `1,2,3,...X,Y,MT`
* POS: mutation position
* ID (optional): mutation ID, default `.`
* REF (optional): reference nucleotide
* ALT (optional): alternative nucleotide
* QUAL (optional): quality score
* FILTER (optional): `PASS` if SNP sites were filtered, default `.`
* INFO (optional): additional information
* FORMAT (optional): recording variant genotype information for the sample
* Sample (optional): the corresponding data in `FORMAT` field.

* example
> \#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample\
> 1	3025531	.	T	A	64.28	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;QD=32.14;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:76,6,0\
> 1	3037125	.	A	C	68.28	PASS	AC=2;AF=1.00;AN=2;DP=2;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;QD=34.14;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:80,6,0\
> 1	5170624	.	A	G	434.6	SnpCluster	AC=1;AF=0.500;AN=2;BaseQRankSum=2.479;DP=17;ExcessHet=3.0103;FS=5.315;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=25.56;ReadPosRankSum=-1.640;SOR=0.662	GT:AD:DP:GQ:PL	0/1:5,12:17:99:442,0,142\
> 1	85864585	.	T	A,C	771.02	PASS	AC=1,1;AF=0.500,0.500;AN=2;DP=20;ExcessHet=3.0103;FS=0.000;MLEAC=1,1;MLEAF=0.500,0.500;MQ=60.00;QD=31.86;SOR=1.022	GT:AD:DP:GQ:PL	1/2:0,5,14:19:99:788,579,564,209,0,167


### 2. BED format file
Contains fields below, more details see [BED format demonstration UCSC](http://genome.ucsc.edu/FAQ/FAQformat#format1)
* \# chr: chromosome number, `1,2,3,...X,Y,MT`
* chromStart: m6A signal start position on chromosome
* chromEnd: m6A signal end position on chromosome
* name: ENSEMBL gene ID
* score: significant score(adjusted p.value), generate by peak calling tools, less than `0.05`
* strand: `+` or `-`
* thickStart: The starting position at which the feature is drawn thickly
* thickEnd: The ending position at which the feature is drawn thickly
* itemRgb: An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. 
* blockCount: sub-block number of the m6A signal peak, integer, `≥1`
* blockSizes: block size of each sub-block, separate by `,`
* blockStarts: block start position on chromosome, separate by `,`

> \# chr	chromStart	chromEnd	name	score	strand	thickStart	thickEnd	itemRgb	blockCount	blockSizes	blockStarts\
> 1	9747647	9747845	ENSMUSG00000097893	7.1e-05	+	9747647	9747845	0	1	198,	0\
> 1	16105773	16105923	ENSMUSG00000025921	4.9e-05	+	16105773	16105923	0	1	150,	0\
> 1	33739519	33739819	ENSMUSG00000004768	0.0032	+	33739519	33739819	0	1	300,	0\
> 1	34180162	34180463	ENSMUSG00000026131	0.00022	+	34180162	34180463	0	1	301,	0\
> 1	34306583	34307612	ENSMUSG00000026131	0.00038	+	34306583	34307612	0	2	68,283,	0,746

### OUTPUT FILE DESCRIPTION
### 1. ASE gene detection output
When the algorithm finishes running, the following files will be in your output folder:
* error.log: error logs generated during program execution
* logout.log: normal logs generated during program execution
* snp_location.txt: the information of SNP loci used in the algorithmic computation process
* aseGene.txt (if you specify an output filename, it will be the name you specified): the result of ASE gene detect
   * geneId: The Ensembl ID of the gene
   * geneName: The name of the gene being tested
   * pValue: The p-value computed by the algorithm
   * qValue: The result of the p-value after undergoing Benjamini-Hochberg (BH) correction
   * snpNum: The number of SNP loci used in the algorithm for identifying the gene
   * major/minorAlleleReads: The number of reads available for calculation on the major/minor allele of the gene
   * majorAlleleFrequency: The frequency of the major allele
   * majorAlleleBase: The major allele base at the SNP loci used for calculation

### 2. ASM m6A signal detection output
When the algorithm finishes running, the following files will be in your output folder:
* error.log: error logs generated during program execution
* logout.log: normal logs generated during program execution
* snp_location.txt: the information of SNP loci used in the ASE algorithmic computation process
* aseRes.txt: the ASE gene identified by the algorithm
* peak_with_snp.txt: the information of SNP loci covered by the m6A peak
* asmPeak.txt (if you specify an output filename, it will be the name you specified): the result of ASM m6A signal detect
   * chr: Chromosome
   * peakStart: The genomic position of the start point of the m6A peak
   * peakEnd: The genomic position of the end point of the m6A peak
   * geneId: The Ensembl ID of the gene to which the m6A peak belongs
   * geneName: The name of the gene to which the m6A peak belongs
   * pValue: The p-value computed by the algorithm
   * qValue: The result of the p-value after undergoing Benjamini-Hochberg (BH) correction
   * snpNum: The number of SNP loci used in the algorithm for identifying the m6A peak
   * major/minorAlleleReads: The number of reads available for calculation on the major/minor allele of the peak
   * majorAlleleFrequency: Major allele frequency, where the major allele is based on the genotype of the INPUT sample. **If this value is less than 0.5**, it indicates allele-specific m6A modification occurring on the minor allele. Conversely, if it's greater than 0.5, it indicates allele-specific m6A modification occurring on the major allele
   * majorAlleleBase: The major allele base at the SNP loci covered by the identified m6A peak

### 3. Sample-Specific ASM m6A signal detection output
When the algorithm finishes running, the following files will be in your output folder:
* error.log: error logs generated during program execution
* logout.log: normal logs generated during program execution
* sample1/snp_location.txt: the information of SNP loci used in the ASE algorithmic computation process for sample1
* sample1/peak_with_snp.txt: the information of SNP loci covered by the m6A peak in sample1
* sample2/snp_location.txt: the information of SNP loci used in the ASE algorithmic computation process for sample2
* sample2/peak_with_snp.txt: the information of SNP loci covered by the m6A peak in sample2
* sampleSpecificAsm6A.txt: the identification results of sample-specific ASM m6A signal
    * chr: Chromosome
    * peakStart: The genomic position of the start point of the m6A peak
    * peakEnd: The genomic position of the end point of the m6A peak
    * geneId: The Ensembl ID of the gene to which the m6A peak belongs
    * geneName: The name of the gene to which the m6A peak belongs
    * sample1MajorFrequency: The major allele frequency of peak m6A in sample 1, where the major allele is based on the genotype of the INPUT sample. **If this value is less than 0.5**, it indicates allele-specific m6A modification occurring on the minor allele. Conversely, if it's greater than 0.5, it indicates allele-specific m6A modification occurring on the major allele. **If the m6A peak has no snp loci in sample1 that can be used for calculation, it will be represented as -.**
    * sample2MajorFrequency: The major allele frequency of peak m6A in sample 2, where the major allele is based on the genotype of the INPUT sample. **If this value is less than 0.5**, it indicates allele-specific m6A modification occurring on the minor allele. Conversely, if it's greater than 0.5, it indicates allele-specific m6A modification occurring on the major allele. **If the m6A peak has no snp loci in sample2 that can be used for calculation, it will be represented as -.**
    * sample1MajorHaplotype: The major allele genotype within the m6A peak in sample1. **If the m6A peak has no snp loci in sample1 that can be used for calculation, it will be represented as -.**
    * sample2MajorHaplotype: The major allele genotype within the m6A peak in sample2. **If the m6A peak has no snp loci in sample2 that can be used for calculation, it will be represented as -.**
    * sample1PValue: 
        * If the m6A peak has no snp loci in sample1 that can be used for calculation, it will be represented as -
        * If the specifically modified haplotypes of the m6A peak are the same in both sample1 and sample2, then the value represents the result recalculated using the data from both samples simultaneously
        * If there is no resampling calculation performed and the peak has SNP loci available for calculation in sample1, then this value is calculated solely based on the data sampled from sample1
    * sample1QValue: The Benjamini-Hochberg corrected value for the sample1PValue
    * sample2PValue: 
        * If the m6A peak has no snp loci in sample2 that can be used for calculation, it will be represented as -
        * If the specifically modified haplotypes of the m6A peak are the same in both sample1 and sample2, then the value represents the result recalculated using the data from both samples simultaneously
        * If there is no resampling calculation performed and the peak has SNP loci available for calculation in sample2, then this value is calculated solely based on the data sampled from sample2
    * sample2QValue: The Benjamini-Hochberg corrected value for the sample2PValue
    * specificSample: In which sample does the peak have sample-specific ASM
        * -: There is no sample-specific m6A modification in either of the two samples
        * sample1: There is sample-specific m6A modification in sample1
        * sample2: There is sample-specific m6A modification in sample2
        * sample1/sample2: There is sample-specific m6A modification in both sample