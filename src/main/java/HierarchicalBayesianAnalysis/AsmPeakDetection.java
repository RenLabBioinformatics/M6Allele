package HierarchicalBayesianAnalysis;

import HeterozygoteSiteAnalysis.HeterozygoteReadsCount;
import HeterozygoteSiteAnalysis.PeakWithSNV;
import Threshold.GPD.GPDFunction;
import Threshold.RecordParamDTO;
import org.apache.commons.cli.*;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;

public class AsmPeakDetection {
    private final String gtfFile;
    private final String asmPeakFile;
    private final String peakCoveredSnpFile;
    private final int ipSNPReadInfimum;
    private final int samplingTime;
    private final int burnIn;
    private final int threadNumber;
    private final Logger log;
    private HashMap<String, HashMap<String, List<String>>> peakSnpReadsCount;
    private final HashMap<String, HashMap<String, Integer>> peakMajorMinorAlleleCount = new HashMap<>();
    private final HashMap<String, HashMap<String, Integer>> peakMajorMinorBackground = new HashMap<>();
    private final HashMap<String, Double> peakMajorAlleleFrequency = new HashMap<>();
    private final Map<String, Double> geneMajorAlleleOddRatio;
    private final HashMap<String, String> peakCoveredGene = new HashMap<>();
    private HashMap<String, List<String[]>> statisticForTest = new HashMap<>();
    private final HashMap<String, LinkedHashSet<String>> peakMajorAlleleNucleotide = new HashMap<>();
    private HashMap<String, HashMap<String, String>> geneNames;
    private final ArrayList<String> asmQValue = new ArrayList<>();
    private ReentrantLock lock;
    public Map<String, RecordParamDTO> distributionType2GpdParam = new HashMap<>();
    public Map<String, List<String>> distributionType2GeneIds = new HashMap<>();
    private final HashMap<String, HashMap<Double, ArrayList<String>>> asmPValue = new HashMap<>();
    private final Map<String, Integer> peakLabel2SnpSize = new HashMap<>();
    private final Map<String, List<String>> peakLabel2ValidSnpInfo = new HashMap<>();
    private final HashMap<String, HashMap<Integer, String[]>> aseGeneIdName2SnpBaseCount;
    
    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param peakBedFile BED format file via MeRIP-seq IP data
     * @param vcfFile VCF format file via MeRIP-seq INPUT data
     * @param asmPeakFile test result output file
     * @param ipSNPReadInfimum reads coverage threshold when filter INPUT sample SNV sites
     * @param samplingTime sampling time, default 10000
     * @param burnIn burn in time, default 2000
     * @param threadNumber thread number, default 2
     * @param log log4j instance
     */
    public AsmPeakDetection(String gtfFile, String peakBedFile, String vcfFile,
                            String asmPeakFile, int ipSNPReadInfimum, int samplingTime, int burnIn,
                            int threadNumber, Logger log, Map<String, Double> geneMajorAlleleOddRatio,
                            List<String> ipBamFile, List<String> ipBaiFile, HashMap<String,
                               HashMap<Integer, String[]>> aseGeneIdName2SnpBaseCount) {
        this.gtfFile = gtfFile;
        this.aseGeneIdName2SnpBaseCount = aseGeneIdName2SnpBaseCount;
        String outputDir = new File(asmPeakFile).getParent();
        this.log = log;
        getGpdFunction();
        getDistributionType2GeneIds();
        this.asmPeakFile = asmPeakFile;
        this.ipSNPReadInfimum = ipSNPReadInfimum;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.threadNumber = threadNumber;
        
        this.log.info("locate SNP record in ip.bam to corresponding m6A signal peaks");
        this.peakCoveredSnpFile = new File(outputDir, "peak_with_snp.txt").getAbsolutePath();
        PeakWithSNV pws = new PeakWithSNV(gtfFile, peakBedFile, vcfFile, this.peakCoveredSnpFile);
        pws.locateSnvInPeak(ipBamFile, ipBaiFile);
        this.log.info("SNP locations in " + this.peakCoveredSnpFile);
        this.geneMajorAlleleOddRatio = geneMajorAlleleOddRatio;
    }
    
    
    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "AsmPeakDetection contains following parameters: ";
        String footer = "";
        
        try {
            commandLine = setCommandLine(args, options);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        }
        
        if (commandLine.hasOption("h")) {
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(0);
        }
        
        // default parameters
        String gtfFile = null, bedFile = null, aseVcfFile = null, outputFile, outputDir;
        int ipSNPCoverageInfimum = 10, samplingTime = 50000, burnIn = 10000, threadNumber = 2;
        List<String> inputBamFilePath = new ArrayList<>();
        List<String> inputBaiFilePath = new ArrayList<>();
        List<String> ipBamFilePath = new ArrayList<>();
        List<String> ipBaiFilePath = new ArrayList<>();
        
        if (!commandLine.hasOption("o")) {
            outputFile = new File(System.getProperty("user.dir"), "asmPeak.txt").getAbsolutePath();
        } else {
            outputFile = commandLine.getOptionValue("o");
        }
        
        outputDir = new File(outputFile).getParent();
        Logger logger = initLog(outputDir);
        
        if (!commandLine.hasOption("bed")) {
            logger.error("Peak calling BED format file can not be empty");
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            File bed = new File(commandLine.getOptionValue("bed"));
            if (!bed.exists() || !bed.isFile()) {
                logger.error("invalid bed file path: " + bed.getAbsolutePath());
                System.exit(2);
            }
            bedFile = bed.getAbsolutePath();
        }
        
        if (!commandLine.hasOption("vcf")) {
            logger.error("SNP calling VCF file can not be empty");
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            File vcf = new File(commandLine.getOptionValue("vcf"));
            if (!vcf.exists() || !vcf.isFile()) {
                logger.error("invalid vcf file path: " + vcf.getAbsolutePath());
                System.exit(2);
            }
            aseVcfFile = vcf.getAbsolutePath();
        }
        File temp;
        if (!commandLine.hasOption("inputBam")) {
            logger.error("inputBamFile can not be empty");
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            String[] inputBamFiles = commandLine.getOptionValue("inputBam").split(",");
            if (inputBamFiles.length == 0) {
                logger.error("doesn't have input bam files!");
                System.exit(2);
            }
            for (String inputBamFile : inputBamFiles) {
                temp = new File(inputBamFile);
                if (!temp.exists() || !temp.isFile()) {
                    logger.error("invalid input bam file path: " + temp.getAbsolutePath());
                    System.exit(2);
                }
                inputBamFilePath.add(temp.getAbsolutePath());
            }
        }
        
        if (!commandLine.hasOption("inputBai")) {
            logger.error("inputBaiFile can not be empty");
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            String[] inputBaiFiles = commandLine.getOptionValue("inputBai").split(",");
            if (inputBaiFiles.length == 0) {
                logger.error("doesn't have input bai files!");
                System.exit(2);
            }
            for (String inputBaiFile : inputBaiFiles) {
                temp = new File(inputBaiFile);
                if (!temp.exists() || !temp.isFile()) {
                    logger.error("invalid input bai file path: " + temp.getAbsolutePath());
                    System.exit(2);
                }
                inputBaiFilePath.add(temp.getAbsolutePath());
            }
        }
        
        if (!commandLine.hasOption("ipBam")) {
            logger.error("ipBamFile can not be empty");
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            String[] ipBamFiles = commandLine.getOptionValue("ipBam").split(",");
            if (ipBamFiles.length == 0) {
                logger.error("doesn't have ip bam files!");
                System.exit(2);
            }
            for (String ipBamFile : ipBamFiles) {
                temp = new File(ipBamFile);
                if (!temp.exists() || !temp.isFile()) {
                    logger.error("invalid Ip bam file path: " + temp.getAbsolutePath());
                    System.exit(2);
                }
                ipBamFilePath.add(temp.getAbsolutePath());
            }
        }
        
        if (!commandLine.hasOption("ipBai")) {
            logger.error("ipBaiFile can not be empty");
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            String[] ipBaiFiles = commandLine.getOptionValue("ipBai").split(",");
            if (ipBaiFiles.length == 0) {
                logger.error("doesn't have ip bai files!");
                System.exit(2);
            }
            for (String ipBaiFile : ipBaiFiles) {
                temp = new File(ipBaiFile);
                if (!temp.exists() || !temp.isFile()) {
                    logger.error("invalid Ip bai file path: " + temp.getAbsolutePath());
                    System.exit(2);
                }
                ipBaiFilePath.add(temp.getAbsolutePath());
            }
        }
        
        if (!commandLine.hasOption("g")) {
            logger.error("GTF format file can not be empty");
            help.printHelp("java -jar M6Allele.jar AsmPeakDetection", header, options, footer, true);
            System.exit(2);
        } else {
            File gtf = new File(commandLine.getOptionValue("g"));
            if (!gtf.exists() || !gtf.isFile()) {
                logger.error("invalid gtf file path: " + gtf.getAbsolutePath());
                System.exit(2);
            }
            gtfFile = gtf.getAbsolutePath();
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
        if (commandLine.hasOption("rc")) {
            ipSNPCoverageInfimum = Integer.parseInt(commandLine.getOptionValue("rc"));
        }
        if (commandLine.hasOption("t")) {
            if (Integer.parseInt(commandLine.getOptionValue("t")) < 0) {
                System.err.println("invalid thread number, should be a positive integer");
                System.exit(2);
            }
            threadNumber = Integer.parseInt(commandLine.getOptionValue("t"));
        }
        
        AseGeneDetection agd = new AseGeneDetection(gtfFile, aseVcfFile, outputFile,
                ipSNPCoverageInfimum, samplingTime, burnIn, threadNumber, logger, inputBamFilePath, inputBaiFilePath);
        agd.getTestResult(0);
        
        AsmPeakDetection apd = new AsmPeakDetection(gtfFile, bedFile, aseVcfFile, outputFile,
                ipSNPCoverageInfimum, samplingTime, burnIn, threadNumber, logger,
                agd.getGeneMajorAlleleOddRatio(), ipBamFilePath, ipBaiFilePath, agd.getGeneAlleleReads());
        apd.getTestResult();
    }
    
    public HashMap<String, String> getPeakCoveredGenes() {
        return this.peakCoveredGene;
    }
    
    public HashMap<String, HashMap<String, String>> getGeneNames() {
        return geneNames;
    }
    
    public HashMap<String, HashMap<String, List<String>>> getPeakSnpReadsCount() {
        return this.peakSnpReadsCount;
    }
    
    public void setPeakSnpReadsCount(HashMap<String, HashMap<String, List<String>>> peakSnpReadsCount) {
        this.peakSnpReadsCount = peakSnpReadsCount;
    }
    
    public void getTestResult() {
        this.parseGTFFile();
        this.asmPeakTest();
        this.bhRecalibrationOfEachPeak();
        this.outputResult();
    }

    /**
     * get major allele and minor allele nucleotide and reads counts of each MeRIP-seq INPUT sample SNV sites covered by m6A peak
     * [chr1:peakStart:peakEnd -> [position1:majorNC, position2:majorNC],....]
     */
    public void getPeakSNPReadsCount() {
        HeterozygoteReadsCount hrc = new HeterozygoteReadsCount(this.peakCoveredSnpFile, this.ipSNPReadInfimum, this.log);
        this.peakSnpReadsCount = hrc.getMajorMinorHaplotype(false);
    }

    /**
     * test the ASM significant p value of a m6A peak
     */
    public void asmPeakTest() {
        // [chr1: [peak1: position1: [major: count, minor: count], position2:[major: count, minor:count]], chr2:....]
        this.getPeakSNPReadsCount();
        this.dataPreparation();
        this.hierarchicalModelTest();
        this.recordValidMajorInfo();
    }

    /**
     * prepare data for hierarchical test
     */
    public void dataPreparation() {
        // pos -> [majorNc:majorCount:bamFileName, minorNc:minorCount:bamFileName....]
        HashMap<String, List<String>> rnaSeqPeakSnvAlleleReads;
        List<String> rnaSeqReads;
        String[] majorCount, minorCount, majorBackground, minorBackground;
        String backgroundMajor, backgroundMinor;
        List<String> rnaSeqMajor, rnaSeqMinor, bgMajor, bgMinor;

        
        // peakLabel   chr:peakStart:peakEnd:geneId
        for (String peakLabel : this.peakSnpReadsCount.keySet()) {
            String chrNum = peakLabel.split(":")[0];
            // SNV sites covered by the m6A signal
            rnaSeqPeakSnvAlleleReads = this.peakSnpReadsCount.get(peakLabel);
            
            // List [majorNc:count:bam1:snp1, majorNc:count:bam1:snp1, majorNc:count:bam2:snp1, majorNc:count:bam2:snp2,...]
            rnaSeqMajor = new ArrayList<>();
            // List [minorNc:count:bam1:snp1, minorNc:count:bam1:snp1, minorNc:count:bam2:snp1, minorNc:count:bam2:snp2,...]
            rnaSeqMinor = new ArrayList<>();
            // List[majorNc:count1:snp1, minor:count2:snp2, ...]
            bgMajor = new ArrayList<>();
            bgMinor = new ArrayList<>();
            for (String position : rnaSeqPeakSnvAlleleReads.keySet()) {
                
                // List[major:count:bam1, minor:count:bam1, major:count:bam2, minor:count:bam2.....]
                rnaSeqReads = rnaSeqPeakSnvAlleleReads.get(position);
                
                for (int i = 0; i < rnaSeqReads.size(); i+=2) {
                    rnaSeqMajor.add(rnaSeqReads.get(i) + ":" + position);
                    rnaSeqMinor.add(rnaSeqReads.get(i + 1) + ":" + position);
                }
                
                backgroundMajor = "no" + ":" + "-1" + ":" + position;
                backgroundMinor = "no" + ":" + "-1" + ":" + position;
                
                bgMajor.add(backgroundMajor);
                bgMinor.add(backgroundMinor);
            }

            if (rnaSeqMajor.size() == 0) {
                continue;
            }
            majorCount = new String[rnaSeqMajor.size()];
            minorCount = new String[rnaSeqMinor.size()];
            majorBackground = new String[bgMajor.size()];
            minorBackground = new String[bgMinor.size()];
            //snpPosition = new int[snpPos.size()];
            assert majorCount.length == minorCount.length;
            assert majorCount.length == majorBackground.length;
            assert majorBackground.length == minorBackground.length;
            for (int i = 0; i < majorCount.length; i++) {
                majorCount[i] = rnaSeqMajor.get(i);
                minorCount[i] = rnaSeqMinor.get(i);
            }
            
            for (int i = 0; i < majorBackground.length; i++) {
                majorBackground[i] = bgMajor.get(i);
                minorBackground[i] = bgMinor.get(i);
            }

            List<String[]> statistic = new ArrayList<>();
            statistic.add(majorCount);
            statistic.add(minorCount);
            statistic.add(majorBackground);
            statistic.add(minorBackground);
            this.statisticForTest.put(peakLabel, statistic);

            int totalMajorBkg, totalMinorBkg;
            totalMajorBkg = 0;
            totalMinorBkg = 0;
            
            HashMap<String, Integer> background = new HashMap<>();
            background.put("major", totalMajorBkg);
            background.put("minor", totalMinorBkg);
            this.peakMajorMinorBackground.put(peakLabel, background);
        }
        
        if (this.statisticForTest.isEmpty()) {
            this.log.error("contains no peaks with SNV sites for hierarchical test, please check the input data");
            System.exit(2);
        }
    }
    
    /**
     * get gpd function to compute threshold
     */
    public void getGpdFunction() {
        InputStream paramInputStream = AseGeneDetection.class.getClassLoader().getResourceAsStream("gpdParameter.txt");
        try (BufferedReader br = new BufferedReader(new InputStreamReader(paramInputStream))){
            String line = br.readLine();
            String[] infos;
            while ((line = br.readLine()) != null) {
                infos = line.split("\t");
                RecordParamDTO dto = new RecordParamDTO();
                dto.gpdGamma = Double.parseDouble(infos[1]);
                dto.gpdSigma = Double.parseDouble(infos[2]);
                dto.N = Integer.parseInt(infos[3]);
                dto.Nt = Integer.parseInt(infos[4]);
                dto.t = Double.parseDouble(infos[5]);
                dto.samplingMafs = Arrays.stream(infos[6].split(",")).mapToDouble(Double::parseDouble).boxed().collect(Collectors.toList());
                distributionType2GpdParam.put(infos[0], dto);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public void getDistributionType2GeneIds() {
        InputStream geneRecords = AsmPeakDetection.class.getClassLoader().getResourceAsStream("distribution2Genes.txt");
        try (BufferedReader br = new BufferedReader(new InputStreamReader(geneRecords))){
            String line;
            String[] infos;
            while ((line = br.readLine()) != null) {
                infos = line.split("\t");
                String distributionType = infos[0];
                List<String> geneIds = Arrays.asList(infos[1].split(","));
                distributionType2GeneIds.put(distributionType, geneIds);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * calculate m6A signal p value via hierarchical model
     */
    public void hierarchicalModelTest() {
        this.log.info("hierarchical Bayesian model test start");
        ExecutorService threadPoolExecutor = Executors.newFixedThreadPool(this.threadNumber);
        this.lock = new ReentrantLock();
        CountDownLatch countDown = new CountDownLatch(this.statisticForTest.size());
        double tenPercent = Math.ceil(this.statisticForTest.size() * 0.1);

        RunTest task = (String name) -> {
            return new Runnable() {
                @Override
                public void run() {
                    List<String[]> statistic;
                    double df, aveDepth, scaleParam;
                    String[] majorCount, minorCount;
                    HierarchicalBayesianModelAsm hb;
                    try {
                        // name ->  chr:peakStart:peakEnd:geneId
                        statistic = statisticForTest.get(name);
                        // List [majorNc:count:bam1:snp1, majorNc:count:bam1:snp1, majorNc:count:bam2:snp1....]
                        majorCount = statistic.get(0);
                        minorCount = statistic.get(1);
                        
                        String[] labelArr = name.split(":");
                        String geneId = labelArr[3];
                        String geneName = geneNames.get(labelArr[0]).getOrDefault(geneId, "unknown");
                        String label = String.join("->", new String[]{geneId, geneName});
                        double ExpressOddRatio;
                        
                        try {
                            ExpressOddRatio = geneMajorAlleleOddRatio.get(label);
                        } catch (Exception e) {
                            return;
                        }
                        // snpPos -> [majorAllele:Count:bamFile1, minorAllele:count:bamFile1,
                        // majorAllele:count:bam2, minorAllele:count:bam2....]
                        HashMap<Integer, String[]> snpPos2BaseCount = aseGeneIdName2SnpBaseCount.get(label);
                        // get p value via hierarchical model
                        if (majorCount.length == 1 && Integer.parseInt(minorCount[0].split(":")[1]) == 0 &&
                                Integer.parseInt(majorCount[0].split(":")[1]) <= 35) {
                            df = 2;
                        } else {
                            df = Math.max(3, majorCount.length);
                        }
                        aveDepth = Arrays.stream(majorCount).map(t -> Double.valueOf(t.split(":")[1])).mapToDouble(t -> t).average().getAsDouble();
                        if (majorCount.length == 1) {
                            scaleParam = 50;
                        } else {
                            scaleParam = (aveDepth - 15 < 0.000001)? 50: 100;
                        }
                        double pVal, peakOddRatio, peakMAF;
                        String distributionType;
                        hb = new HierarchicalBayesianModelAsm(df, scaleParam, samplingTime,
                                burnIn, majorCount, minorCount, ExpressOddRatio, snpPos2BaseCount);
                        boolean flag = hb.testSignificant();
                        
                        if (!flag) {
                            return;
                        }
                        peakOddRatio = Math.exp(hb.quantifyGeneLOR());
                        peakMAF = Math.min(1.0, peakOddRatio / (peakOddRatio + 1));
                        PValueParam param = getPValue(peakMAF, geneId);
                        pVal = param.pvalue;
                        distributionType = param.disType;
                        lock.lock();
                        
                        peakLabel2SnpSize.put(name, hb.getSnpNum());
                        peakLabel2ValidSnpInfo.put(name, hb.getValidSnpInfo());
                        HashMap<Double, ArrayList<String>> pvalue2PeakLabel = asmPValue.getOrDefault(distributionType, new HashMap<>());
                        ArrayList<String> samePValPeaks = pvalue2PeakLabel.getOrDefault(pVal, new ArrayList<>());
                        samePValPeaks.add(name);
                        pvalue2PeakLabel.put(pVal, samePValPeaks);
                        asmPValue.put(distributionType, pvalue2PeakLabel);
                        peakMajorAlleleFrequency.put(name, peakMAF);
                    } catch (Exception e) {
                        log.error("error occurs on record " + name);
                        log.error(e.getMessage());
                    } finally {
                        countDown.countDown();
                        if (countDown.getCount() % tenPercent == 0) {
                            double proportion = 100 - 10.0 * countDown.getCount() / tenPercent;
                            if (proportion >= 0)
                                log.info(proportion + "% completed");
                        }
                        if (lock.isHeldByCurrentThread()) {
                            lock.unlock();
                        }
                    }
                }
            };
        };

        this.log.debug(this.statisticForTest.size() + " m6A signal peaks to be tested");
        for (String peakLabel : this.statisticForTest.keySet()) {
            Runnable runnable = task.runTask(peakLabel);
            threadPoolExecutor.execute(runnable);
        }
        try {
            countDown.await();
        } catch (InterruptedException ie) {
            this.log.error("analysis interrupted");
            this.log.error(ie.getMessage());
        } finally {
            this.statisticForTest = null;
            try {
                threadPoolExecutor.shutdown();
                if (!threadPoolExecutor.awaitTermination(1000, TimeUnit.MILLISECONDS))
                    threadPoolExecutor.shutdownNow();
            } catch (InterruptedException ie) {
                threadPoolExecutor.shutdownNow();
            }
            this.lock = null;
        }
        this.log.info("model asm test complete");
    }
    
    private String getDistributionType(String geneId) {
        String distributionType = null;
        for(Map.Entry<String, List<String>> entry : distributionType2GeneIds.entrySet()) {
            if (entry.getValue().contains(geneId)) {
                distributionType = entry.getKey();
                break;
            }
        }
        return distributionType;
    }
    
    public PValueParam getPValue(double predictAsm, String geneId) {
        predictAsm = Math.max(predictAsm, 1 - predictAsm);
        String distributionType = getDistributionType(geneId);
        double pvalue;
        
        if (distributionType == null) {
            distributionType = "unknown";
        }
        RecordParamDTO paramDTO = distributionType2GpdParam.getOrDefault(distributionType, null);
        GPDFunction gpdFunction = new GPDFunction(paramDTO.gpdGamma, paramDTO.gpdSigma);
        if (predictAsm >= paramDTO.t) {
            pvalue = ((double) paramDTO.Nt / paramDTO.N) * (1 - gpdFunction.cdf(predictAsm - paramDTO.t));
        } else {
            double count = 0;
            for (int i = 0; i < paramDTO.samplingMafs.size(); i++) {
                if (paramDTO.samplingMafs.get(i) >= predictAsm) {
                    count++;
                }
            }
            pvalue = count / paramDTO.samplingMafs.size();
        }
        return new PValueParam(pvalue, distributionType);
    }

    /**
     * recalibrate the p value with BH method, get significant q value
     */
    public void bhRecalibrationOfEachPeak() {
        this.log.info("start recalibrating p values of hierarchical model");
        this.log.info("sorting test result in order");
        for (String distributionType : asmPValue.keySet()) {
            HashMap<Double, ArrayList<String>> pValue2PeakLabel = this.asmPValue.get(distributionType);
            
            ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(pValue2PeakLabel.entrySet());
    
            // sort p value from large to small
            sortedPVals.sort((o1, o2) -> o2.getKey().compareTo(o1.getKey()));
            
            int totalPeak = 0;
            for (Map.Entry<Double, ArrayList<String>> sortedPVal : sortedPVals) {
                totalPeak += sortedPVal.getValue().size();
            }
    
            int rankage = totalPeak;
            double prevQValue = 1.0, qValue;
            String pValString, qValString;
            // entry   pvalue -> [peakLabel1, peakLabel2....]
            for (Map.Entry<Double, ArrayList<String>> entry: sortedPVals) {
                Double pVal = entry.getKey();
                ArrayList<String> samePValPeaks = entry.getValue();
                // sort items with its SNV number when p value is same
                HashMap<String, Integer> samePValPeaksSNVs = new HashMap<>();
                for (String peak: samePValPeaks) {
                    samePValPeaksSNVs.put(peak, this.peakLabel2SnpSize.get(peak));
                }
        
                // sort items with its major allele frequency when p value and SNV numbers are same
                HashMap<String, Double> samePValPeakMajorAlleleFrequency = new HashMap<>();
                for (String peak: samePValPeaks) {
                    samePValPeakMajorAlleleFrequency.put(peak, this.peakMajorAlleleFrequency.get(peak));
                }
        
                List<Map.Entry<String, Integer>> samePValPeakEntry = new ArrayList<>(samePValPeaksSNVs.entrySet());
                samePValPeakEntry.sort((o1, o2) -> {
            
                    String peak1 = o1.getKey(), peak2 = o2.getKey();
                    Integer peak1SNVs = o1.getValue(), peak2SNVs = o2.getValue();
                    Double peak1MAF = samePValPeakMajorAlleleFrequency.get(peak1), peak2MAF = samePValPeakMajorAlleleFrequency.get(peak2);
                    if (peak1SNVs.equals(peak2SNVs)) {
                        return peak2MAF.compareTo(peak1MAF);
                    } else {
                        return peak2SNVs.compareTo(peak1SNVs);
                    }
                });
        
                for (Map.Entry<String, Integer> geneEntry: samePValPeakEntry) {
                    String peak = geneEntry.getKey();
                    qValue = Math.min(prevQValue, pVal * totalPeak / rankage);
                    if (qValue - prevQValue < 0.00001) {
                        prevQValue = qValue;
                    }
                    rankage--;
            
                    pValString = String.valueOf(pVal);
                    qValString = String.valueOf(qValue);
                    this.asmQValue.add(String.join("->", new String[]{peak, pValString, qValString}));
                }
            }
        }
        this.log.info("recalibration complete.");
    }
    
    public void recordValidMajorInfo() {
        // peakLabel2ValidSnpInfo  chr:peakStart:peakEnd:geneId -> [majorNc:cnt:bam1:snpPos, minorNc:cnt:bam1:snpPos, ...]
        for (String label : peakLabel2ValidSnpInfo.keySet()) {
            int majorNum = 0;
            int minorNum = 0;
            List<String> snpInfos = peakLabel2ValidSnpInfo.get(label);
            LinkedHashSet<String> snpPosMajorNc = new LinkedHashSet<>();
            for (int i = 0; i < snpInfos.size(); i += 2) {
                String[] majorInfos = snpInfos.get(i).split(":");
                String[] minorInfos = snpInfos.get(i + 1).split(":");
                majorNum += Integer.parseInt(majorInfos[1]);
                minorNum += Integer.parseInt(minorInfos[1]);
                snpPosMajorNc.add(String.join(":", new String[]{majorInfos[3], majorInfos[0]}));
            }
            HashMap<String, Integer> record = new HashMap<>();
            record.put("major", majorNum);
            record.put("minor", minorNum);
            peakMajorMinorAlleleCount.put(label, record);
            peakMajorAlleleNucleotide.put(label, snpPosMajorNc);
        }
    }

    /**
     * write the test result into file
     */
    private void outputResult() {
        ArrayList<String[]> outputRecord = new ArrayList<>();
        BufferedWriter bfw = null;
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(this.asmPeakFile)));
            String line, label, geneId, geneName, chrNum, peakStart, peakEnd, majorAlleleReads, minorAlleleReads,
                    majorAlleleBackground, minorAlleleBackground, pValue, qValue;
            LinkedHashSet<String> majorAlleleNucleotides;
            String[] info, rec, finalInfo;
            String majorAlleleFrequency;
            int snvNum;
            bfw.write("chr\tpeakStart\tpeakEnd\tgeneId\tgeneName\tpValue\tqValue\tsnpNum\tmajor/minorAlleleReads\tmajorAlleleFrequency\tmajorAlleleBase\n");
            for (String record: this.asmQValue) {
                rec = record.split("->");
                label = rec[0];
                // info = [chrNum, peakStart, peakEnd]
                info = label.split(":");
                chrNum = info[0];
                peakStart = info[1];
                peakEnd = info[2];
                geneId = info[3];
                pValue = rec[1];
                qValue = rec[2];
                //geneId = this.peakCoveredGene.get(label);
                geneName = this.geneNames.get(chrNum).getOrDefault(geneId, "unknown");
                majorAlleleFrequency = String.valueOf(this.peakMajorAlleleFrequency.get(label));
                majorAlleleNucleotides = this.peakMajorAlleleNucleotide.get(label);
                String majorAlleleRecords = this.getString(majorAlleleNucleotides);
                //snvNum = this.peakSNVNum.get(label);
                snvNum = majorAlleleNucleotides.size();
                majorAlleleReads = String.valueOf(this.peakMajorMinorAlleleCount.get(label).get("major"));
                minorAlleleReads = String.valueOf(this.peakMajorMinorAlleleCount.get(label).get("minor"));
                majorAlleleBackground = String.valueOf(this.peakMajorMinorBackground.get(label).get("major"));
                minorAlleleBackground = String.valueOf(this.peakMajorMinorBackground.get(label).get("minor"));
                finalInfo = new String[]{chrNum, peakStart, peakEnd, geneId, geneName, pValue, qValue,
                        Integer.toString(snvNum), majorAlleleReads + "," + minorAlleleReads,
                        majorAlleleBackground+","+minorAlleleBackground, majorAlleleFrequency,
                        majorAlleleRecords};
                outputRecord.add(finalInfo);
            }

            outputRecord.sort((o1, o2) -> {
                Double q1 = Double.parseDouble(o1[6]), q2 = Double.parseDouble(o2[6]);
                if (!q1.equals(q2))
                    return q1.compareTo(q2);
                // sort records with same q-value after BH recalibration by SNV number
                Integer snvCount1 = Integer.parseInt(o1[7]), snvCount2 = Integer.parseInt(o2[7]);
                return snvCount2.compareTo(snvCount1);
            });
            for (String[] record: outputRecord) {
                line = String.join("\t", record[0], record[1], record[2], record[3], record[4], record[5], record[6],
                        record[7], record[8], record[10], record[11]);
                bfw.write(line);
                bfw.newLine();
            }
            this.log.info("result file " + this.asmPeakFile);
        } catch (IOException ie) {
            this.log.error(ie.getMessage());
        } finally {
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }
    
    public Map<String, String[]> mergeTestRes() {
        Map<String, String[]> peakLabel2OutRes = new HashMap<>();
        String label, geneId, geneName, chrNum, peakStart, peakEnd, majorAlleleReads, minorAlleleReads,
                pValue, qValue;
        LinkedHashSet<String> majorAlleleNucleotides;
        String[] info, rec, finalInfo;
        String majorAlleleFrequency;
        int snvNum;
        // record  chr:peakStart:peakEnd:geneId->pvalue->qvalue
        for (String record: this.asmQValue) {
            rec = record.split("->");
            label = rec[0];
            // info = [chrNum, peakStart, peakEnd, geneId]
            info = label.split(":");
            chrNum = info[0];
            peakStart = info[1];
            peakEnd = info[2];
            geneId = info[3];
            pValue = rec[1];
            qValue = rec[2];
            geneName = this.geneNames.get(chrNum).getOrDefault(geneId, "unknown");
            majorAlleleFrequency = String.valueOf(this.peakMajorAlleleFrequency.get(label));
            majorAlleleNucleotides = this.peakMajorAlleleNucleotide.get(label);
            String majorAlleleRecords = this.getString(majorAlleleNucleotides);
            snvNum = majorAlleleNucleotides.size();
            majorAlleleReads = String.valueOf(this.peakMajorMinorAlleleCount.get(label).get("major"));
            minorAlleleReads = String.valueOf(this.peakMajorMinorAlleleCount.get(label).get("minor"));
            finalInfo = new String[]{chrNum, peakStart, peakEnd, geneId, geneName, pValue, qValue,
                    Integer.toString(snvNum), majorAlleleReads + "," + minorAlleleReads,
                    majorAlleleFrequency, majorAlleleRecords};
            peakLabel2OutRes.put(label, finalInfo);
        }
        return peakLabel2OutRes;
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

    public void parseGTFFile() {
        BufferedReader bfr = null;
        this.geneNames = new HashMap<>();
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.gtfFile)));
            String line = "", chrNum, geneId, geneName;
            String[] info, geneInfo;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    if (!info[2].equalsIgnoreCase("gene")) {
                        continue;
                    }
                    chrNum = info[0];
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];
                    HashMap<String, String> chrGenes = this.geneNames.getOrDefault(chrNum, new HashMap<>());
                    chrGenes.put(geneId, geneName);
                    this.geneNames.put(chrNum, chrGenes);
                }
            }
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(2);
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private String[] getGeneInfo(String recordInfo) {
        String[] info = recordInfo.split("; ");
        String geneName = "unknown", geneId = null;
        for (String s: info) {
            if (s.startsWith("gene_id")) {
                String[] name = s.split(" ");
                geneId = name[1].substring(1, name[1].length() -1);
            }
            if (s.startsWith("gene_name")) {
                String[] name = s.split(" ");
                geneName = name[1].substring(1, name[1].length() -1);
            }
        }

        return new String[] {geneId, geneName};
    }

    /**
     * initial log4j Logger instance
     * @param logHome output directory of log file
     * @return Logger instance
     */
    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(AsmPeakDetection.class);
    }

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("vcf", "vcf_file", true, "VCF format file generate by RNA-seq or MeRIP-seq data SNP calling process, required");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("bed", "peak_bed_file", true, "Peak calling output result in BED format, required");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("inputBam", "input_bam_file", true, "input bam files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("inputBai", "input_bai_file", true, "input bam file index files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("ipBam", "ip_bam_file", true, "ip bam files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("ipBai", "ip_bai_file", true, "ip bam file index files, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("g", "gtf", true, "GTF format file, required");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("o", "output", true, "ASM m6A signal test output file. Optional, default ./asmPeak.txt");
        option.setRequired(false);
        options.addOption(option);

        option = new Option("rc", "reads_coverage", true, "reads coverage threshold using for filter RNA-seq or MeRIP-seq data SNVs in VCF file (aim for reducing FDR). Optional, default 5");
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

        option = new Option("h", "help", false, "help message of AsmPeakDetection");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
