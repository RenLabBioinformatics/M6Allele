package DifferentiationAnalysis;

import HierarchicalBayesianAnalysis.AseGeneDetection;
import HierarchicalBayesianAnalysis.AsmPeakDetection;
import org.apache.commons.cli.*;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.log4j.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author: tanglin
 * @date: 2023/03/17 10:58
 * @Description:
 */
public class PairedSampleASM {
    private final Logger logger;
    private final AseGeneDetection sample1AseGeneDetector;
    private final AseGeneDetection sample2AseGeneDetector;
    private final String mergeBedFile;
    private final String sample1VcfFile;
    private final String sample2VcfFile;
    private final String sample1ResFile;
    private final String sample2ResFile;
    private final List<String> sample1InputBamFiles;
    private final List<String> sample1IpBamFiles;
    private final List<String> sample2IpBamFiles;
    private final List<String> sample1IpBaiFiles;
    private final List<String> sample2IpBaiFiles;
    private final String resFile;
    
    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param sample1VcfFile sample1 VCF format file via MeRIP-seq INPUT data
     * @param sample2VcfFile sample2 VCF format file via MeRIP-seq INPUT data
     * @param outputDir result output directory
     * @param readsCoverageThreshold reads coverage threshold when filter INPUT sample SNV sites, default 10
     * @param samplingTime sampling time, default 50000
     * @param burnIn burn in time, default 10000
     * @param threadNumber thread number, default 2
     * @param logger log4j instance
     */
    public PairedSampleASM(String gtfFile, String sample1VcfFile, String sample2VcfFile,
                           String mergePeakBedFile, String outputDir, int readsCoverageThreshold,
                           int samplingTime, int burnIn, int threadNumber, Logger logger, List<String> sample1InputBamFiles,
                           List<String> sample1InputBaiFiles, List<String> sample2InputBamFiles, List<String> sample2InputBaiFiles,
                           List<String> sample1IpBamFiles, List<String> sample1IpBaiFiles, List<String> sample2IpBamFiles,
                           List<String> sample2IpBaiFiles) {
        this.logger = logger;
        this.sample1InputBamFiles = sample1InputBamFiles;
        this.sample1IpBamFiles = sample1IpBamFiles;
        this.sample1IpBaiFiles = sample1IpBaiFiles;
        this.sample2IpBamFiles = sample2IpBamFiles;
        this.sample2IpBaiFiles = sample2IpBaiFiles;
        this.mergeBedFile = mergePeakBedFile;
        this.sample1VcfFile = sample1VcfFile;
        this.sample2VcfFile = sample2VcfFile;
        File sample1OutputDir = new File(outputDir, "sample1");
        if (!sample1OutputDir.exists()) {
            sample1OutputDir.mkdirs();
        }
        File sample2OutputDir = new File(outputDir, "sample2");
        if (!sample2OutputDir.exists()) {
            sample2OutputDir.mkdirs();
        }
        this.sample1ResFile = new File(sample1OutputDir, "sample1Asm6ARes.txt").getAbsolutePath();
        this.sample1AseGeneDetector = new AseGeneDetection(gtfFile, sample1VcfFile,
                sample1ResFile, readsCoverageThreshold,
                samplingTime, burnIn, threadNumber, logger, sample1InputBamFiles, sample1InputBaiFiles);
        
        this.sample2ResFile = new File(sample2OutputDir, "sample2Asm6ARes.txt").getAbsolutePath();
        this.sample2AseGeneDetector = new AseGeneDetection(gtfFile, sample2VcfFile,
                sample2ResFile, readsCoverageThreshold,
                samplingTime, burnIn, threadNumber, logger, sample2InputBamFiles,
                sample2InputBaiFiles);
        resFile = new File(outputDir, "sampleSpecificAsm6A.txt").getAbsolutePath();
    }
    
    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "SampleSpecificASM contains following parameters: ";
        String footer = "";
        
        try {
            commandLine = setCommandLine(args, options);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar M6Allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        }
        
        if (commandLine.hasOption("h")) {
            help.printHelp("java -jar M6Allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(0);
        }
        
        // default parameters
        String gtfFile = null, sample1VcfFile = null, sample2VcfFile = null, outputDir, mergeBedFile = null;
        int samplingTime = 50000, burnIn = 10000, readsCoverageThreshold = 10, threadNumber = 2;
        double threshold = 0.05;
        List<String> sample1InputBamFilePath = new ArrayList<>();
        List<String> sample1InputBaiFilePath = new ArrayList<>();
        List<String> sample2InputBamFilePath = new ArrayList<>();
        List<String> sample2InputBaiFilePath = new ArrayList<>();
    
        List<String> sample1IpBamFilePath = new ArrayList<>();
        List<String> sample1IpBaiFilePath = new ArrayList<>();
        List<String> sample2IpBamFilePath = new ArrayList<>();
        List<String> sample2IpBaiFilePath = new ArrayList<>();
        
        if (commandLine.hasOption("threshold")) {
            threshold = Double.parseDouble(commandLine.getOptionValue("threshold"));
        }
        if (!commandLine.hasOption("o")) {
            outputDir = new File(System.getProperty("user.dir")).getAbsolutePath();
        }
        else {
            outputDir = commandLine.getOptionValue("o");
        }
        Logger logger = initLog(outputDir);
        
        if (!commandLine.hasOption("g")) {
            logger.error("GTF annotation file can not be empty");
            help.printHelp("java -jar M6Allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File gtf = new File(commandLine.getOptionValue("g"));
            if (!gtf.exists() || !gtf.isFile()) {
                logger.error("invalid file path: " + gtf.getAbsolutePath());
                System.exit(2);
            }
            gtfFile = gtf.getAbsolutePath();
        }
        
        if (!commandLine.hasOption("s1Vcf")) {
            logger.error("sample1 SNP calling VCF file can not be empty");
            help.printHelp("java -jar M6Allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("s1Vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid s1Vcf file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            sample1VcfFile = vcf.getAbsolutePath();
        }
        
        if (!commandLine.hasOption("s2Vcf")) {
            logger.error("sample2 SNP calling VCF file can not be empty");
            help.printHelp("java -jar M6Allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("s2Vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid s2Vcf file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            sample2VcfFile = vcf.getAbsolutePath();
        }
    
        if (!commandLine.hasOption("bed")) {
            logger.error("Peak calling BED format file can not be empty");
            help.printHelp("java -jar M6Allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            File bed = new File(commandLine.getOptionValue("bed"));
            if (!bed.exists() || !bed.isFile()) {
                logger.error("invalid s1Bed file path: " + bed.getAbsolutePath());
                System.exit(2);
            }
            mergeBedFile = bed.getAbsolutePath();
        }
        
        
        if (commandLine.hasOption("rc")) {
            readsCoverageThreshold = Integer.parseInt(commandLine.getOptionValue("rc"));
        }
        if (commandLine.hasOption("s")) {
            samplingTime = Integer.parseInt(commandLine.getOptionValue("s"));
        }
        if (commandLine.hasOption("b")) {
            burnIn = Integer.parseInt(commandLine.getOptionValue("b"));
        }
        if (samplingTime <= 500 || burnIn <= 100) {
            logger.error("sampling times larger than 500 and burn in times at least 100");
            System.exit(2);
        }
        if (commandLine.hasOption("t")) {
            if (Integer.parseInt(commandLine.getOptionValue("t")) <= 0) {
                System.err.println("invalid thread number, should be a positive integer");
                System.exit(2);
            }
            threadNumber = Integer.parseInt(commandLine.getOptionValue("t"));
        }
        
        readBamAndBaiFiles(commandLine, help, header, options, footer, "s1InputBam", "s1InputBai", sample1InputBamFilePath, sample1InputBaiFilePath, logger);
        readBamAndBaiFiles(commandLine, help, header, options, footer, "s2InputBam", "s2InputBai", sample2InputBamFilePath, sample2InputBaiFilePath, logger);
        readBamAndBaiFiles(commandLine, help, header, options, footer, "s1IpBam", "s1IpBai", sample1IpBamFilePath, sample1IpBaiFilePath, logger);
        readBamAndBaiFiles(commandLine, help, header, options, footer, "s2IpBam", "s2IpBai", sample2IpBamFilePath, sample2IpBaiFilePath, logger);
        
        PairedSampleASM ssasm = new PairedSampleASM(gtfFile, sample1VcfFile, sample2VcfFile,
                mergeBedFile, outputDir, readsCoverageThreshold, samplingTime,
                burnIn, threadNumber, logger, sample1InputBamFilePath, sample1InputBaiFilePath, sample2InputBamFilePath, sample2InputBaiFilePath,
                sample1IpBamFilePath, sample1IpBaiFilePath, sample2IpBamFilePath, sample2IpBaiFilePath);
        ssasm.testPairedSampleAsm(gtfFile, readsCoverageThreshold, samplingTime, burnIn, threadNumber, threshold);
    }
    
    private static void readBamAndBaiFiles(CommandLine commandLine, HelpFormatter help, String header, Options options, String footer,
                                           String bam, String bai, List<String> bamFilePaths, List<String> baiFilePaths, Logger logger) {
        File temp;
        if (!commandLine.hasOption(bam)) {
            logger.error("bam file can not be empty");
            help.printHelp("java -jar M6Allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            String[] bamFiles = commandLine.getOptionValue(bam).split(",");
            if (bamFiles.length == 0) {
                logger.error(bam + " doesn't have bam files!");
                System.exit(2);
            }
            for (String bamFile : bamFiles) {
                temp = new File(bamFile);
                if (!temp.exists() || !temp.isFile()) {
                    logger.error("invalid bam file path: " + temp.getAbsolutePath());
                    System.exit(2);
                }
                bamFilePaths.add(temp.getAbsolutePath());
            }
        }
        
        if (!commandLine.hasOption(bai)) {
            logger.error("bai file can not be empty");
            help.printHelp("java -jar M6Allele.jar SampleSpecificASM", header, options, footer, true);
            System.exit(2);
        } else {
            String[] baiFiles = commandLine.getOptionValue(bai).split(",");
            if (baiFiles.length == 0) {
                logger.error(bai + " doesn't have bai files!");
                System.exit(2);
            }
            for (String baiFile : baiFiles) {
                temp = new File(baiFile);
                if (!temp.exists() || !temp.isFile()) {
                    logger.error("invalid bai file path: " + temp.getAbsolutePath());
                    System.exit(2);
                }
                baiFilePaths.add(temp.getAbsolutePath());
            }
        }
        
        if (bamFilePaths.size() != baiFilePaths.size()) {
            logger.error("the number of bam file and the number of bai file are not equal");
            System.exit(2);
        }
    }
    
    public void testPairedSampleAsm(String gtfFile, int readsCoverageThreshold,
                                    int samplingTime, int burnIn, int threadNumber, double threshold) {
        HashMap<String, HashMap<Integer, String[]>> s1GeneAlleleReads = this.sample1AseGeneDetector.getGeneAlleleReads();
        HashMap<String, HashMap<Integer, String[]>> s2GeneAlleleReads = this.sample2AseGeneDetector.getGeneAlleleReads();
        
        Set<String> s1Genes = new HashSet<>(s1GeneAlleleReads.keySet());
        Set<String> s2Genes = new HashSet<>(s2GeneAlleleReads.keySet());
        
        Set<String> commonGeneIdNames = new HashSet<>(s1GeneAlleleReads.keySet());
        commonGeneIdNames.retainAll(s2Genes);
        Set<String> commonGeneIds = commonGeneIdNames.stream().map(geneIdName -> geneIdName.split("->")[0])
                .collect(Collectors.toSet());
    
        Set<String> onlyS1GeneIdNames = new HashSet<>(s1GeneAlleleReads.keySet());
        onlyS1GeneIdNames.removeAll(s2Genes);
        Set<String> onlyS1GeneIds = onlyS1GeneIdNames.stream().map(geneIdName -> geneIdName.split("->")[0])
                .collect(Collectors.toSet());
        
        s2Genes.removeAll(s1Genes);
        Set<String> onlyS2GeneIds = s2Genes.stream().map(geneIdName -> geneIdName.split("->")[0])
                .collect(Collectors.toSet());
        
        
        this.geneHaplotype();
        this.logger.info("compute sample1 ase background...");
        this.sample1AseGeneDetector.dataPreparation();
        this.sample1AseGeneDetector.aseGeneTest();
        this.logger.info("compute sample2 ase background...");
        this.sample2AseGeneDetector.dataPreparation();
        this.sample2AseGeneDetector.aseGeneTest();
        this.logger.info("sample ase background compute ending...");
        
        this.logger.info("get sample1 and sample2's peak reads count...");
        AsmPeakDetection sample1AsmPeakDetection = new AsmPeakDetection(gtfFile, mergeBedFile, this.sample1VcfFile,
                this.sample1ResFile, readsCoverageThreshold, samplingTime, burnIn, threadNumber,
                logger, this.sample1AseGeneDetector.getGeneMajorAlleleOddRatio(), sample1IpBamFiles, sample1IpBaiFiles,
                this.sample1AseGeneDetector.getGeneAlleleReads());
        sample1AsmPeakDetection.parseGTFFile();
        sample1AsmPeakDetection.getPeakSNPReadsCount();
        AsmPeakDetection sample2AsmPeakDetection = new AsmPeakDetection(gtfFile, mergeBedFile, sample2VcfFile,
                sample2ResFile, readsCoverageThreshold, samplingTime, burnIn, threadNumber,
                logger, sample2AseGeneDetector.getGeneMajorAlleleOddRatio(), sample2IpBamFiles, sample2IpBaiFiles,
                sample2AseGeneDetector.getGeneAlleleReads());
        sample2AsmPeakDetection.parseGTFFile();
        sample2AsmPeakDetection.getPeakSNPReadsCount();
        
        
        this.peakHaplotype(sample1AsmPeakDetection, sample2AsmPeakDetection);
        this.logger.info("compute sample1 asm...");
        // chr:peakStart:peakEnd:geneId -> {chrNum, peakStart, peakEnd, geneId, geneName, pValue, qValue,
        //                    Integer.toString(snvNum), majorAlleleReads + "," + minorAlleleReads,
        //                    majorAlleleFrequency, majorAlleleHaplotype}
        Map<String, String[]> sample1Output = this.getTestResult(sample1AsmPeakDetection);
        this.logger.info("sample1 asm compute success...");
        this.logger.info("compute sample2 asm...");
        Map<String, String[]> sample2Output = this.getTestResult(sample2AsmPeakDetection);
        this.logger.info("sample1 asm compute success...");
        
        Map<String, List<String>> validOutput = new HashMap<>();
        
        Set<String> commonPeakLabel = new HashSet<>();
        
        this.processSingleSampleOutputRes(sample1Output, commonPeakLabel, onlyS1GeneIds, validOutput, 0, threshold);
        
        this.processSingleSampleOutputRes(sample2Output, commonPeakLabel, onlyS2GeneIds, validOutput, 1, threshold);
        
        commonPeakLabel.removeIf(e -> !commonGeneIds.contains(e.split(":")[3]));
        
        Set<String> needResamplingPeakLabel = new HashSet<>();
        for (String label : commonPeakLabel) {
            String[] sample1Infos = sample1Output.get(label);
            String[] sample2Infos = sample2Output.get(label);
            if (sample1Infos == null || sample2Infos == null) {
                continue;
            }
            double sample1QValue = Double.parseDouble(sample1Infos[6]);
            double sample2QValue = Double.parseDouble(sample2Infos[6]);
            // chrNum, peakStart, peakEnd, geneId, geneName sample1MajorFrequency sample2MajorFrequency
            // sample1MajorHaplotype sample2MajorHaplotype p1Value q1Value p2Value q2Value specificSample
            List<String> validInfo = new ArrayList<>();
            validInfo.add(sample1Infos[0]);
            validInfo.add(sample1Infos[1]);
            validInfo.add(sample1Infos[2]);
            validInfo.add(sample1Infos[3]);
            validInfo.add(sample1Infos[4]);
            validInfo.add(sample1Infos[9]);
            validInfo.add(sample2Infos[9]);
            validInfo.add(sample1Infos[10]);
            validInfo.add(sample2Infos[10]);
            validInfo.add(sample1Infos[5]);
            validInfo.add(sample1Infos[6]);
            validInfo.add(sample2Infos[5]);
            validInfo.add(sample2Infos[6]);
            if (sample1QValue >= threshold && sample2QValue >= threshold) {
                validInfo.add("-");
            } else if (sample1QValue < threshold && sample2QValue >= threshold) {
                validInfo.add("sample1");
            } else if (sample1QValue >= threshold && sample2QValue < threshold) {
                validInfo.add("sample2");
            } else if (!isEqualsHaplotype(sample1Infos[10], sample2Infos[10])){
                validInfo.add("sample1/sample2");
            } else {
                needResamplingPeakLabel.add(label);
            }
            validOutput.put(label, validInfo);
        }
        
        
        SpecificAsmDetection specificAsmDetection = new SpecificAsmDetection(sample1AsmPeakDetection.getPeakSnpReadsCount(),
                sample2AsmPeakDetection.getPeakSnpReadsCount(), sample1AseGeneDetector.getGeneMajorAlleleOddRatio(),
                sample2AseGeneDetector.getGeneMajorAlleleOddRatio(), sample1AseGeneDetector.getGeneAlleleReads(),
                sample2AseGeneDetector.getGeneAlleleReads(), needResamplingPeakLabel, logger, threadNumber,
                sample1AsmPeakDetection.getGeneNames(), samplingTime, burnIn);
        
        //this.peakHaplotype(sample1AsmPeakDetection, sample2AsmPeakDetection, needResamplingPeakLabel);
        specificAsmDetection.asmPairedCompare();
        // label-> pVal -> qVal
        List<String> pairedCompareQValue = specificAsmDetection.getPairedCompareQValue();
        Map<String, LinkedHashSet<String>> peakMajorNc = specificAsmDetection.getPeakMajorNc();
        //Map<String, Double> specificAsmLogOddRatio = specificAsmDetection.getPeakMajorLogOddRatio();
        for (String qValueInfos : pairedCompareQValue) {
            String[] infos = qValueInfos.split("->");
            List<String> validInfo = validOutput.get(infos[0]);
            LinkedHashSet<String> majorNc = peakMajorNc.get(infos[0]);
            String majorNcStr = this.getString(majorNc);
            
            validInfo.set(validInfo.size() - 1, infos[2]);
            validInfo.set(validInfo.size() - 2, infos[1]);
            validInfo.set(validInfo.size() - 3, infos[2]);
            validInfo.set(validInfo.size() - 4, infos[1]);
            validInfo.set(validInfo.size() - 5, majorNcStr);
            validInfo.set(validInfo.size() - 6, majorNcStr);
            
            if (Double.parseDouble(infos[2]) >= threshold) {
                validInfo.add("-");
            } else {
                double s1Maf = Double.parseDouble(validInfo.get(5));
                double s2Maf = Double.parseDouble(validInfo.get(6));
                if (s1Maf < 0.5 && s2Maf < 0.5) {
                    if (s1Maf < s2Maf) {
                        validInfo.add("sample1");
                    } else {
                        validInfo.add("sample2");
                    }
                } else if (s1Maf > 0.5 && s2Maf > 0.5) {
                    if (s1Maf < s2Maf) {
                        validInfo.add("sample2");
                    } else {
                        validInfo.add("sample1");
                    }
                } else {
                    validInfo.add("sample1/sample2");
                }
            }
            validOutput.put(infos[0], validInfo);
        }
        
        List<String> discardPeaks = specificAsmDetection.getDiscardPeak();
        for (String discardPeak : discardPeaks) {
            validOutput.remove(discardPeak);
        }
        this.logger.info("analysis complete");
        outputResult(validOutput);
    }
    
    private String getString(LinkedHashSet<String> list) {
        List<String> sortedList = list.stream().sorted((o1, o2) -> {
            Integer pos1 = Integer.parseInt(o1.split(":")[0]);
            Integer pos2 = Integer.parseInt(o2.split(":")[0]);
            return pos1.compareTo(pos2);
        }).collect(Collectors.toList());
        
        String[] str = new String[sortedList.size()];
        for (int i = 0; i < sortedList.size(); i++) {
            str[i] = sortedList.get(i);
        }
        
        return String.join(";", str);
    }
    
    private boolean isEqualsHaplotype(String s1MajorHaplotype, String s2MajorHaplotype) {
        Set<String> s1Infos = Arrays.stream(s1MajorHaplotype.split(";")).collect(Collectors.toSet());
        Set<String> s2Infos = Arrays.stream(s2MajorHaplotype.split(";")).collect(Collectors.toSet());
        Set<String> commonTypes = new HashSet<>(s1Infos);
        commonTypes.retainAll(s2Infos);
        return commonTypes.size() >= Math.min(s1Infos.size(), s2Infos.size()) * 3 / 4;
    }
    
    // status  0: Process the sample1 result  1: Process the sample2 result
    private void processSingleSampleOutputRes(Map<String, String[]> sampleOutputRes, Set<String> commonPeakLabel,
                                              Set<String> geneIdsOnlyInSingSample, Map<String, List<String>> validOutput,
                                              int status, double threshold) {
        for (String peakLabel : sampleOutputRes.keySet()) {
            String[] infos = peakLabel.split(":");
            if (!geneIdsOnlyInSingSample.contains(infos[3])) {
                commonPeakLabel.add(peakLabel);
                continue;
            }
            // chr:peakStart:peakEnd:geneId -> {chrNum, peakStart, peakEnd, geneId, geneName, pValue, qValue,
            //                    Integer.toString(snvNum), majorAlleleReads + "," + minorAlleleReads,
            //                    majorAlleleFrequency, majorAlleleHaplotype}
            String[] onlySingleSampleInfos = sampleOutputRes.get(peakLabel);
            if (Double.parseDouble(onlySingleSampleInfos[6]) >= threshold) {
                continue;
            }
            // chrNum, peakStart, peakEnd, geneId, geneName sample1MajorFrequency sample2MajorFrequency
            // sample1MajorHaplotype sample2MajorHaplotype p1Value q1Value p2Value q2Value specificSample
            List<String> validInfo = new ArrayList<>();
            // 构造有效信息 sample2相关的都用-占位符
            validInfo.add(onlySingleSampleInfos[0]);
            validInfo.add(onlySingleSampleInfos[1]);
            validInfo.add(onlySingleSampleInfos[2]);
            validInfo.add(onlySingleSampleInfos[3]);
            validInfo.add(onlySingleSampleInfos[4]);
            if (status == 0) {
                validInfo.add(onlySingleSampleInfos[9]);
                validInfo.add("-");
                validInfo.add(onlySingleSampleInfos[10]);
                validInfo.add("-");
                validInfo.add(onlySingleSampleInfos[5]);
                validInfo.add(onlySingleSampleInfos[6]);
                validInfo.add("-");
                validInfo.add("-");
                validInfo.add("sample1");
            } else {
                validInfo.add("-");
                validInfo.add(onlySingleSampleInfos[9]);
                validInfo.add("-");
                validInfo.add(onlySingleSampleInfos[10]);
                validInfo.add("-");
                validInfo.add("-");
                validInfo.add(onlySingleSampleInfos[5]);
                validInfo.add(onlySingleSampleInfos[6]);
                validInfo.add("sample2");
            }
            validOutput.put(peakLabel, validInfo);
        }
    }
    
    private void outputResult(Map<String, List<String>> validOutput) {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(resFile))){
            // chrNum, peakStart, peakEnd, geneId, geneName sample1MajorFrequency sample2MajorFrequency
            // sample1MajorHaplotype sample2MajorHaplotype pValue qValue specificSample
            bw.write("chr\tpeakStart\tpeakEnd\tgeneId\tgeneName\tsample1MajorFrequency\tsample2MajorFrequency\t" +
                    "sample1MajorHaplotype\tsample2MajorHaplotype\tsample1PValue\tsample1QValue\tsample2PValue\tsample2QValue\tspecificSample");
            bw.newLine();
            for (String label : validOutput.keySet()) {
                List<String> infos = validOutput.get(label);
                bw.write(String.join("\t", infos));
                bw.newLine();
            }
        } catch (IOException e) {
            this.logger.error(e);
        }
    }
    
    private Map<String, String[]> getTestResult(AsmPeakDetection asmPeakDetectionV2) {
        asmPeakDetectionV2.dataPreparation();
        asmPeakDetectionV2.hierarchicalModelTest();
        asmPeakDetectionV2.recordValidMajorInfo();
        asmPeakDetectionV2.bhRecalibrationOfEachPeak();
        return asmPeakDetectionV2.mergeTestRes();
    }
    
    private void geneHaplotype() {
        int bamFileNum = this.sample1InputBamFiles.size();
        // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count:bam1, minorAllele:count:bam1], pos2: ...}, ...}
        HashMap<String, HashMap<Integer, String[]>> sample1GeneAlleleReads = this.sample1AseGeneDetector.getGeneAlleleReads();
        HashMap<String, HashMap<Integer, String[]>> sample2GeneAlleleReads = this.sample2AseGeneDetector.getGeneAlleleReads();
        
        HashSet<String> sample1Genes = new HashSet<>(sample1GeneAlleleReads.keySet());
        HashSet<String> sample2Genes = new HashSet<>(sample2GeneAlleleReads.keySet());
        sample1Genes.retainAll(sample2Genes);
        if (sample1Genes.isEmpty()) {
            this.logger.info("don't detect common genes between two samples");
            return;
        } else {
            this.logger.info("start unifying haplotype information");
        }
        /*sample1GeneAlleleReads.keySet().removeIf(key -> !sample1Genes.contains(key));
        sample2GeneAlleleReads.keySet().removeIf(key -> !sample1Genes.contains(key));*/
        
        
        for (String gene : sample1Genes) {
            // pos1: [majorAllele:count:bam1, minorAllele:count:bam1], pos2: ...}, ...
            HashMap<Integer, String[]> sample1Snp2BaseInfo = sample1GeneAlleleReads.get(gene);
            HashMap<Integer, String[]> sample2Snp2BaseInfo = sample2GeneAlleleReads.get(gene);
            Set<Integer> sample1SnpPoss = new HashSet<>(sample1Snp2BaseInfo.keySet());
            Set<Integer> sample2SnpPoss = new HashSet<>(sample2Snp2BaseInfo.keySet());
            sample1SnpPoss.retainAll(sample2SnpPoss);
            sample1Snp2BaseInfo.keySet().removeIf(key -> !sample1SnpPoss.contains(key));
            sample2Snp2BaseInfo.keySet().removeIf(key -> !sample1SnpPoss.contains(key));
            
            for (Integer snpPos : sample1SnpPoss) {
                // majorNc:count:bam1, minor:count:bam1, majorNc:count:bam2, minorNc:count:bam2...
                String[] sample1OriginBaseInfos = sample1Snp2BaseInfo.get(snpPos);
                String[] sample2OriginBaseInfos = sample2Snp2BaseInfo.get(snpPos);
                List<String> sample1NewBaseInfoList = new ArrayList<>();
                List<String> sample2NewBaseInfoList = new ArrayList<>();
                for (int i = 0; i < bamFileNum; i++) {
                    String bamMark = "bam" + (i + 1);
                    // originSampleBaseInfo : [Major:count:bam1}, Minor:count:bam1]
                    List<String> originSample1BaseInfo = getSampleOriginBaseInfoInBam(sample1OriginBaseInfos, bamMark);
                    List<String> originSample2BaseInfo = getSampleOriginBaseInfoInBam(sample2OriginBaseInfos, bamMark);
                    
                    if (CollectionUtils.isEmpty(originSample1BaseInfo) || CollectionUtils.isEmpty(originSample2BaseInfo)) {
                        continue;
                    }
                    String sample1MajorNc = originSample1BaseInfo.get(0).split(":")[0];
                    String sample1MinorNc = originSample1BaseInfo.get(1).split(":")[0];
                    String sample2MajorNc = originSample2BaseInfo.get(0).split(":")[0];
                    String sample2MinorNc = originSample2BaseInfo.get(1).split(":")[0];
                    if (sample1MajorNc.equalsIgnoreCase(sample2MajorNc) && sample1MinorNc.equalsIgnoreCase(sample2MinorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.addAll(originSample2BaseInfo);
                    } else if (sample1MajorNc.equalsIgnoreCase(sample2MinorNc) && sample1MinorNc.equalsIgnoreCase(sample2MajorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(1));
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(0));
                    } else if (sample1MajorNc.equalsIgnoreCase(sample2MajorNc) && !sample1MinorNc.equalsIgnoreCase(sample2MinorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(0));
                        sample2NewBaseInfoList.add(setCountTo0(originSample1BaseInfo.get(1)));
                    } else if (sample1MajorNc.equalsIgnoreCase(sample2MinorNc) && !sample1MinorNc.equalsIgnoreCase(sample2MajorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(1));
                        sample2NewBaseInfoList.add(setCountTo0(originSample1BaseInfo.get(1)));
                    } else if (!sample1MajorNc.equalsIgnoreCase(sample2MajorNc) && sample1MinorNc.equalsIgnoreCase(sample2MinorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(setCountTo0(originSample1BaseInfo.get(0)));
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(1));
                    } else if (!sample1MajorNc.equalsIgnoreCase(sample2MinorNc) && sample1MinorNc.equalsIgnoreCase(sample2MajorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(setCountTo0(originSample1BaseInfo.get(0)));
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(0));
                    }
                }
                if (sample1NewBaseInfoList.size() != 0) {
                    sample1Snp2BaseInfo.put(snpPos, sample1NewBaseInfoList.toArray(new String[0]));
                } else {
                    sample1Snp2BaseInfo.remove(snpPos);
                }
                if (sample2NewBaseInfoList.size() != 0) {
                    sample2Snp2BaseInfo.put(snpPos, sample2NewBaseInfoList.toArray(new String[0]));
                } else {
                    sample2Snp2BaseInfo.remove(snpPos);
                }
            }
            if (sample1Snp2BaseInfo.isEmpty()) {
                sample1GeneAlleleReads.remove(gene);
            } else {
                sample1GeneAlleleReads.put(gene, sample1Snp2BaseInfo);
            }
            if (sample2Snp2BaseInfo.isEmpty()) {
                sample2GeneAlleleReads.remove(gene);
            } else {
                sample2GeneAlleleReads.put(gene, sample2Snp2BaseInfo);
            }
        }
        this.sample1AseGeneDetector.setGeneAlleleReads(sample1GeneAlleleReads);
        this.sample2AseGeneDetector.setGeneAlleleReads(sample2GeneAlleleReads);
        this.logger.info("unify haplotype completed");
    }
    
    
    private void peakHaplotype(AsmPeakDetection sample1AsmPeakDetection, AsmPeakDetection sample2AsmPeakDetection) {
        // chr:peakStart:peakEnd:geneId -> {pos1: [majorNc:cnt:bam1, minorNc:cnt:bam1, majorNc:cnt:bam2, minorNc:cnt:bam2...], pos2:...}
        HashMap<String, HashMap<String, List<String>>> sample1PeakSnpReadsCount = sample1AsmPeakDetection.getPeakSnpReadsCount();
        HashMap<String, HashMap<String, List<String>>> sample2PeakSnpReadsCount = sample2AsmPeakDetection.getPeakSnpReadsCount();
        int bamFileNum = this.sample1IpBamFiles.size();
        Set<String> sample1Peaks = new HashSet<>(sample1PeakSnpReadsCount.keySet());
        Set<String> sample2Peaks = new HashSet<>(sample2PeakSnpReadsCount.keySet());
        sample1Peaks.retainAll(sample2Peaks);
        
        for (String peak : sample1Peaks) {
            // pos1: [majorAllele:count:bam1, minorAllele:count:bam1], pos2: ...}, ...
            HashMap<String, List<String>> sample1Snp2BaseInfo = sample1PeakSnpReadsCount.get(peak);
            HashMap<String, List<String>> sample2Snp2BaseInfo = sample2PeakSnpReadsCount.get(peak);
            Set<String> sample1SnpPoss = new HashSet<>(sample1Snp2BaseInfo.keySet());
            Set<String> sample2SnpPoss = new HashSet<>(sample2Snp2BaseInfo.keySet());
            sample1SnpPoss.retainAll(sample2SnpPoss);
            sample1Snp2BaseInfo.keySet().removeIf(key -> !sample1SnpPoss.contains(key));
            sample2Snp2BaseInfo.keySet().removeIf(key -> !sample1SnpPoss.contains(key));
            
            for (String snpPos : sample1SnpPoss) {
                // majorNc:count:bam1, minor:count:bam1, majorNc:count:bam2, minorNc:count:bam2...
                List<String> sample1OriginBaseInfos = sample1Snp2BaseInfo.get(snpPos);
                List<String> sample2OriginBaseInfos = sample2Snp2BaseInfo.get(snpPos);
                List<String> sample1NewBaseInfoList = new ArrayList<>();
                List<String> sample2NewBaseInfoList = new ArrayList<>();
                for (int i = 0; i < bamFileNum; i++) {
                    String bamMark = "bam" + (i + 1);
                    // originSampleBaseInfo : [Major:count:bam1, Minor:count:bam1]
                    List<String> originSample1BaseInfo = getSampleOriginBaseInfoInBam(sample1OriginBaseInfos.toArray(new String[0]), bamMark);
                    List<String> originSample2BaseInfo = getSampleOriginBaseInfoInBam(sample2OriginBaseInfos.toArray(new String[0]), bamMark);
                    if (CollectionUtils.isEmpty(originSample1BaseInfo) || CollectionUtils.isEmpty(originSample2BaseInfo)) {
                        continue;
                    }
                    String sample1MajorNc = originSample1BaseInfo.get(0).split(":")[0];
                    String sample1MinorNc = originSample1BaseInfo.get(1).split(":")[0];
                    String sample2MajorNc = originSample2BaseInfo.get(0).split(":")[0];
                    String sample2MinorNc = originSample2BaseInfo.get(1).split(":")[0];
                    if (sample1MajorNc.equalsIgnoreCase(sample2MajorNc) && sample1MinorNc.equalsIgnoreCase(sample2MinorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.addAll(originSample2BaseInfo);
                    } else if (sample1MajorNc.equalsIgnoreCase(sample2MinorNc) && sample1MinorNc.equalsIgnoreCase(sample2MajorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(1));
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(0));
                    } else if (sample1MajorNc.equalsIgnoreCase(sample2MajorNc) && !sample1MinorNc.equalsIgnoreCase(sample2MinorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(0));
                        sample2NewBaseInfoList.add(setCountTo0(originSample1BaseInfo.get(1)));
                    } else if (sample1MajorNc.equalsIgnoreCase(sample2MinorNc) && !sample1MinorNc.equalsIgnoreCase(sample2MajorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(1));
                        sample2NewBaseInfoList.add(setCountTo0(originSample1BaseInfo.get(1)));
                    } else if (!sample1MajorNc.equalsIgnoreCase(sample2MajorNc) && sample1MinorNc.equalsIgnoreCase(sample2MinorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(setCountTo0(originSample1BaseInfo.get(0)));
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(1));
                    } else if (!sample1MajorNc.equalsIgnoreCase(sample2MinorNc) && sample1MinorNc.equalsIgnoreCase(sample2MajorNc)) {
                        sample1NewBaseInfoList.addAll(originSample1BaseInfo);
                        sample2NewBaseInfoList.add(setCountTo0(originSample1BaseInfo.get(0)));
                        sample2NewBaseInfoList.add(originSample2BaseInfo.get(0));
                    }
                }
                if (sample1NewBaseInfoList.size() != 0) {
                    sample1Snp2BaseInfo.put(snpPos, sample1NewBaseInfoList);
                } else {
                    sample1Snp2BaseInfo.remove(snpPos);
                }
                if (sample2NewBaseInfoList.size() != 0) {
                    sample2Snp2BaseInfo.put(snpPos, sample2NewBaseInfoList);
                } else {
                    sample2Snp2BaseInfo.remove(snpPos);
                }
            }
            if (sample1Snp2BaseInfo.isEmpty()) {
                sample1PeakSnpReadsCount.remove(peak);
            } else {
                sample1PeakSnpReadsCount.put(peak, sample1Snp2BaseInfo);
            }
            if (sample2Snp2BaseInfo.isEmpty()) {
                sample2PeakSnpReadsCount.remove(peak);
            } else {
                sample2PeakSnpReadsCount.put(peak, sample2Snp2BaseInfo);
            }
        }
        sample1AsmPeakDetection.setPeakSnpReadsCount(sample1PeakSnpReadsCount);
        sample2AsmPeakDetection.setPeakSnpReadsCount(sample2PeakSnpReadsCount);
        this.logger.info("unifying peak haplotype information end...");
        
    }
    
    private List<String> getSampleOriginBaseInfoInBam(String[] sampleOriginBaseInfos, String bamMark) {
        // return [MajorInfoDTO, MinorInfoDTO]
        List<String> res = new ArrayList<>();
        for (String sampleOriginBaseInfo : sampleOriginBaseInfos) {
            // sampleOriginBaseInfo: majorAllele:count:bam1
            if (sampleOriginBaseInfo.contains(bamMark)) {
                res.add(sampleOriginBaseInfo);
            }
        }
        return res;
    }
    
    private String setCountTo0(String originSnpInfo) {
        String[] tempMajorInfo = originSnpInfo.split(":");
        tempMajorInfo[1] = "0";
        return String.join(":", tempMajorInfo);
    }
    
    
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(PairedSampleASM.class);
    }
    
    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("s1Vcf", "sample1VcfFile", true, "VCF format file generate by sample1 RNA-seq or MeRIP-seq data SNP calling process, required");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("s2Vcf", "sample2VcfFile", true, "VCF format file generate by sample2 RNA-seq or MeRIP-seq data SNP calling process, required");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("bed", "mergePeakBedFile", true, "common bed format peak file generate by sample1 and sample2 MeRIP-seq data");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("g", "gtf", true, "GTF annotation file, required");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("o", "output", true, "ASE gene test output file. Optional, default ./");
        option.setRequired(false);
        options.addOption(option);
        
        option = new Option("s1InputBam", "s1InputBamFiles", true, "sample1 input bam files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("s1InputBai", "s1InputBaiFiles", true, "sample1 input bam file index files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("s2InputBam", "s2InputBamFiles", true, "sample2 input bam files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("s2InputBai", "s2InputBaiFiles", true, "sample2 input bam file index files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("s1IpBam", "s1IpBamFiles", true, "sample1 ip bam files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("s1IpBai", "s1IpBaiFiles", true, "sample1 ip bam file index files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("s2IpBam", "s2IpBamFiles", true, "sample2 ip bam files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("s2IpBai", "s2IpBaiFiles", true, "sample2 ip bam file index files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("rc", "reads_coverage", true, "reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR). Optional, default 10");
        option.setRequired(false);
        options.addOption(option);
        
        option = new Option("s", "sampling", true, "sampling times, larger than 500. Optional, default 50000");
        option.setRequired(false);
        options.addOption(option);
        
        option = new Option("b", "burn", true, "burn-in times, more than 100 and less than sampling times. Optional, default 10000");
        option.setRequired(false);
        options.addOption(option);
        
        option = new Option("t", "thread", true, "thread number for running test. Optional, default 2");
        option.setRequired(false);
        options.addOption(option);
        
        option = new Option("h", "help", false, "help message of SampleSpecificASM");
        option.setRequired(false);
        options.addOption(option);
    
        option = new Option("threshold", "significantThreshold", false, "the threshold for determining whether it is a sample-specific Asm6A modification, default 0.05");
        option.setRequired(false);
        options.addOption(option);
        
        CommandLineParser parser = new DefaultParser();
        
        return parser.parse(options, args);
    }
    

}
