# M6Allele 1.0

## HARDWARE/SOFTWARE REQUIREMENTS
* Java 1.8
* Windows / Linux / Mac OS

## INSTALLATION
* clone the repo,
```
git clone https://github.com/Jakob666/allele-specificM6A.git
```
* target JAR package
```
cd ./allele-specificM6A
```
make sure the directory contains `M6Allele.jar`.

## USAGE
### Overview
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

## FORMAT DECLARATION
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

## OUTPUT FILE DESCRIPTION
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