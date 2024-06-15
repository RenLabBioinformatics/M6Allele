package HierarchicalBayesianAnalysis;

import HeterozygoteSiteAnalysis.GeneSNVRecord;
import HeterozygoteSiteAnalysis.VcfSnpMatchGene;
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

/**
 * Test ASE genes with VCF format file
 */
public class AseGeneDetection {
    private String aseGeneFile;
    // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count:bam1, minorAllele:count:bam1],...}, ...}
    private HashMap<String, HashMap<Integer, String[]>> geneAlleleReads;
    private final HashMap<String, HashMap<Integer, String[]>> geneBackgroundReads;
    private final HashMap<String, HashMap<String, Integer>> geneReadsCount = new HashMap<>();
    private final HashMap<String, HashMap<String, Integer>> geneBackgroundCount = new HashMap<>();
    private final int samplingTime;
    private final int burnIn;
    private final int threadNumber;
    //private HashMap<Double, ArrayList<String>> geneAsePValue = new HashMap<>();
    private final HashMap<String, HashMap<Double, ArrayList<String>>> geneAsePValue = new HashMap<>();
    private final HashMap<String, Integer> genesSNVs = new HashMap<>();
    private final HashMap<String, Double> geneMajorAlleleFrequency = new HashMap<>();
    private final Map<String, Double> geneMajorAlleleOddRatio = new HashMap<>();
    private final HashMap<String, String> geneMajorNucleotide = new HashMap<>();
    private HashMap<String, ArrayList<int[]>> statisticForTest = new HashMap<>();
    private ArrayList<String> geneAseQValue;
    private ReentrantLock lock;
    private final Logger logger;
    public Map<String, RecordParamDTO> distributionType2GpdParam = new HashMap<>();
    public Map<String, List<String>> distributionType2GeneIds = new HashMap<>();

    /**
     * Constructor
     * @param gtfFile GTF annotation file
     * @param vcfFile VCF format file via MeRIP-seq INPUT data
     * @param aseGeneFile test result output file
     * @param readsCoverageThreshold reads coverage threshold when filter INPUT sample SNV sites
     * @param samplingTime sampling time, default 50000
     * @param burnIn burn in time, default 10000
     * @param threadNumber thread number, default 2
     * @param logger log4j instance
     */
    public AseGeneDetection(String gtfFile, String vcfFile, String aseGeneFile,
                            int readsCoverageThreshold, int samplingTime, int burnIn, int threadNumber,
                            Logger logger, List<String> bamFile, List<String> baiFile) {
        
        getGpdFunction();
        getDistributionType2GeneIds();
        VcfSnpMatchGene vsmg = new VcfSnpMatchGene(vcfFile, gtfFile, readsCoverageThreshold);
        
        this.logger = logger;
        this.logger.info("start processing bam file");
        vsmg.parseBamFile(bamFile, baiFile);
        this.geneAlleleReads = vsmg.getGeneAlleleReads();
        this.geneBackgroundReads = null;
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.aseGeneFile = aseGeneFile;
        if (this.aseGeneFile == null) {
            this.logger.error("output file path can not be null");
            System.exit(2);
        }
        this.logger.info("locate SNP record in input.bam to corresponding genes");
        String snvLocationFile = new File(new File(aseGeneFile).getParent(), "snp_location.txt").getAbsolutePath();
        GeneSNVRecord gsr = new GeneSNVRecord(gtfFile, vcfFile, snvLocationFile);
        gsr.locateSnv(bamFile, baiFile);
        this.threadNumber = threadNumber;
        this.logger.info("SNP locations in " + snvLocationFile);
    }

    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        HelpFormatter help = new HelpFormatter();
        String header = "AseGeneDetection contains following parameters: ";
        String footer = "";
        try {
            commandLine = setCommandLine(args, options);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            help.printHelp("java -jar M6Allele.jar AseGeneDetection", header, options, footer, true);
            System.exit(2);
        }

        if (commandLine.hasOption("h")) {
            help.printHelp("java -jar M6Allele.jar AseGeneDetection", header, options, footer, true);
            System.exit(0);
        }

        // default parameters
        String gtfFile = null, aseVcfFile = null, outputFile, outputDir;
        int samplingTime = 50000, burn_in = 10000, readsCoverageThreshold = 10, threadNumber = 2;
        List<String> bamFilePath = new ArrayList<>();
        List<String> baiFilePath = new ArrayList<>();

        if (!commandLine.hasOption("o")) {
            outputFile = new File(System.getProperty("user.dir"), "aseGene.txt").getAbsolutePath();
        } else {
            outputFile = commandLine.getOptionValue("o");
        }
        outputDir = new File(outputFile).getParent();
        Logger logger = initLog(outputDir);

        if (!commandLine.hasOption("g")) {
            logger.error("GTF annotation file can not be empty");
            help.printHelp("java -jar M6Allele.jar AseGeneDetection", header, options, footer, true);
            System.exit(2);
        } else {
            File gtf = new File(commandLine.getOptionValue("g"));
            if (!gtf.exists() || !gtf.isFile()) {
                logger.error("invalid gtf file path: " + gtf.getAbsolutePath());
                System.exit(2);
            }
            gtfFile = gtf.getAbsolutePath();
        }

        if (!commandLine.hasOption("vcf")) {
            logger.error("SNP calling VCF file can not be empty");
            help.printHelp("java -jar M6Allele.jar AseGeneDetection", header, options, footer, true);
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
        if (!commandLine.hasOption("bam")) {
            logger.error("bam file can not be empty");
            help.printHelp("java -jar M6Allele.jar AseGeneDetection", header, options, footer, true);
            System.exit(2);
        } else {
            String[] bamFiles = commandLine.getOptionValue("bam").split(",");
            if (bamFiles.length == 0) {
                logger.error("doesn't have bam files!");
                System.exit(2);
            }
            for (String bamFile : bamFiles) {
                temp = new File(bamFile);
                if (!temp.exists() || !temp.isFile()) {
                    logger.error("invalid bam file path: " + temp.getAbsolutePath());
                    System.exit(2);
                }
                bamFilePath.add(temp.getAbsolutePath());
            }
        }
    
        if (!commandLine.hasOption("bai")) {
            logger.error("bai file can not be empty");
            help.printHelp("java -jar M6Allele.jar AseGeneDetection", header, options, footer, true);
            System.exit(2);
        } else {
            String[] baiFiles = commandLine.getOptionValue("bai").split(",");
            if (baiFiles.length == 0) {
                logger.error("doesn't have bai files!");
                System.exit(2);
            }
            for (String baiFile : baiFiles) {
                temp = new File(baiFile);
                if (!temp.exists() || !temp.isFile()) {
                    logger.error("invalid bai file path: " + temp.getAbsolutePath());
                    System.exit(2);
                }
                baiFilePath.add(temp.getAbsolutePath());
            }
        }
        
        if (bamFilePath.size() != baiFilePath.size()) {
            logger.error("the number of bam file and the number of bai file are not equal");
            System.exit(2);
        }
        

        if (commandLine.hasOption("rc")) {
            readsCoverageThreshold = Integer.parseInt(commandLine.getOptionValue("rc"));
        }
        if (commandLine.hasOption("s")) {
            samplingTime = Integer.parseInt(commandLine.getOptionValue("s"));
        }
        if (commandLine.hasOption("b")) {
            burn_in = Integer.parseInt(commandLine.getOptionValue("b"));
        }
        if (samplingTime <= 500 || burn_in <= 100) {
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
        
        AseGeneDetection agd = new AseGeneDetection(gtfFile, aseVcfFile, outputFile,
                readsCoverageThreshold, samplingTime, burn_in, threadNumber, logger, bamFilePath, baiFilePath);
        agd.getTestResult(1);
    }

    public HashMap<String, String> getGeneMajorNucleotide() {
        return geneMajorNucleotide;
    }

    public HashMap<String, HashMap<Integer, String[]>> getGeneAlleleReads() {
        return this.geneAlleleReads;
    }
    
    public void setGeneAlleleReads(HashMap<String, HashMap<Integer, String[]>> geneAlleleReads) {
        this.geneAlleleReads = geneAlleleReads;
    }
    
    // statusCode: 1 Called by the ase detection algorithm, need to output the ase result
    // otherwise, called by the asm detection algorithm,
    // it needs to output the result of ase but give the output file name itself
    public void getTestResult(int statusCode) {
        this.dataPreparation();
        this.aseGeneTest();
        this.bhRecalibrationOfEachGene();
        this.outputResult(statusCode);
        
    }

    /**
     * data preparation for hierarchical test
     */
    public void dataPreparation() {
        HashMap<Integer, String[]> geneSNVs;
        ArrayList<Integer> majorAlleleCount, minorAlleleCount, majorBackgroundCount, minorBackgroundCount;
        int[] majorCount, minorCount, majorBackground, minorBackground;
        String geneId;
        // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count, minorAllele: count]}, ...}
        this.logger.info("calculate LOR for genes to be tested");
        for (String label : this.geneAlleleReads.keySet()) {

            geneId = label.split("->")[0];
            majorAlleleCount = new ArrayList<>();
            minorAlleleCount = new ArrayList<>();
            majorBackgroundCount = new ArrayList<>();
            minorBackgroundCount = new ArrayList<>();

            // {pos1: [majorAllele:count:bam1, minorAllele:count:bam1,....], pos2: [majorAllele:count:bam1, minorAllele:count:bam1,...]...}
            geneSNVs = this.geneAlleleReads.get(label);
            int snpNum = 0;
            // record the genome positions of major alleles and the corresponding nucleotides
            LinkedList<String> majorNcRecords = new LinkedList<>();
            for (Integer mutPosition : geneSNVs.keySet()) {
                boolean flag = false;
                // [majorAllele:count, minorAllele: count1]
                String[] nucleotideReadsCount = geneSNVs.get(mutPosition);
                String majorAlleleRecord, minorAlleleRecord, majorNC;
                for (int i = 0; i < nucleotideReadsCount.length; i+=2) {
                    majorAlleleRecord = nucleotideReadsCount[i];
                    minorAlleleRecord = nucleotideReadsCount[i + 1];
                    majorNC = majorAlleleRecord.split(":")[0];
                    int major = Integer.parseInt(majorAlleleRecord.split(":")[1]);
                    int minor = Integer.parseInt(minorAlleleRecord.split(":")[1]);
                    majorAlleleCount.add(major);
                    minorAlleleCount.add(minor);
                    if (!flag) {
                        majorNcRecords.add(String.join(":", new String[]{Integer.toString(mutPosition), majorNC}));
                        flag = true;
                    }
                    int alleleCount = (major + minor) / 2;
                    majorBackgroundCount.add(alleleCount);
                    minorBackgroundCount.add(alleleCount);
                }
                snpNum ++;
            }
            if (majorAlleleCount.size() == 0) {
                continue;
            }
            this.genesSNVs.put(label, snpNum);
            this.geneMajorNucleotide.put(geneId, this.getString(majorNcRecords));

            HashMap<String, Integer> readsCount = new HashMap<>();
            readsCount.put("major", this.getSum(majorAlleleCount));
            readsCount.put("minor", this.getSum(minorAlleleCount));
            this.geneReadsCount.put(label, readsCount);
            HashMap<String, Integer> backgroundCount = new HashMap<>();
            backgroundCount.put("major", this.getSum(majorBackgroundCount));
            backgroundCount.put("minor", this.getSum(minorBackgroundCount));
            this.geneBackgroundCount.put(label, backgroundCount);

            assert majorAlleleCount.size() == minorAlleleCount.size();
            assert majorAlleleCount.size() == majorBackgroundCount.size();
            majorCount = new int[majorAlleleCount.size()];
            minorCount = new int[minorAlleleCount.size()];
            majorBackground = new int[majorBackgroundCount.size()];
            minorBackground = new int[minorBackgroundCount.size()];
            for (int i = 0; i < majorAlleleCount.size(); i++) {
                majorCount[i] = majorAlleleCount.get(i);
                minorCount[i] = minorAlleleCount.get(i);

                majorBackground[i] = majorBackgroundCount.get(i);
                minorBackground[i] = minorBackgroundCount.get(i);
            }

            ArrayList<int[]> statistic = new ArrayList<>(4);
            statistic.add(majorCount);
            statistic.add(minorCount);
            statistic.add(majorBackground);
            statistic.add(minorBackground);
            this.statisticForTest.put(label, statistic);
        }
        if (this.statisticForTest.isEmpty()) {
            this.logger.error("contains no genes with SNV sites for hierarchical test, please check the input data");
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
        InputStream geneRecords = AseGeneDetection.class.getClassLoader().getResourceAsStream("distribution2Genes.txt");
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
     * get ASE significant p value with hierarchical model
     */
    public void aseGeneTest() {
        // test ASE gene with Hierarchical model
        this.logger.info("hierarchical Bayesian model test start");
        // init thread pool and lock
        ExecutorService threadPoolExecutor = Executors.newFixedThreadPool(this.threadNumber);
        this.lock = new ReentrantLock();
        CountDownLatch countDown = new CountDownLatch(this.statisticForTest.size());
        long tenPercent = Math.round(this.statisticForTest.size() * 0.1);

        RunTest task = (String name) -> {
            return new Runnable() {
                @Override
                public void run() {
                    ArrayList<int[]> statistic;
                    double df = 5, scaleParam = 10;
                    int[] majorCount, minorCount, majorBackground, minorBackground;
                    HierarchicalBayesianModelAse hb;
                    try {
                        // name  geneId->geneName
                        statistic = statisticForTest.get(name);
                        majorCount = statistic.get(0);
                        minorCount = statistic.get(1);
                        majorBackground = statistic.get(2);
                        minorBackground = statistic.get(3);
                        hb = new HierarchicalBayesianModelAse(df, scaleParam, samplingTime, burnIn,
                                majorCount, minorCount, majorBackground, minorBackground);
                        double p, geneOddRatio, geneMAF;
                        String distributionType;
                        hb.testSignificant();
                        geneOddRatio = Math.exp(hb.quantifyGeneLOR());
                        geneMAF = Math.min(1.0, geneOddRatio / (1 + geneOddRatio));
                        PValueParam param = getPValue(geneMAF, name.split("->")[0]);
                        p = param.pvalue;
                        distributionType = param.disType;
                        if (p == -1) {
                            return;
                        }
                        lock.lock();
                        HashMap<Double, ArrayList<String>> pvalue2GeneIdName = geneAsePValue.getOrDefault(distributionType, new HashMap<>());
                        ArrayList<String> samePValGenes = pvalue2GeneIdName.getOrDefault(p, new ArrayList<>());
                        samePValGenes.add(name);
                        pvalue2GeneIdName.put(p, samePValGenes);
                        geneAsePValue.put(distributionType, pvalue2GeneIdName);
                        geneMajorAlleleFrequency.put(name, geneMAF);
                        geneMajorAlleleOddRatio.put(name, geneOddRatio);
                    } catch (Exception e) {
                        logger.error("error occurs on record " + name);
                        logger.error(e.getMessage());
                    } finally {
                        countDown.countDown();
                        if (countDown.getCount() % tenPercent == 0) {
                            double proportion = 100 - 10.0 * countDown.getCount() / tenPercent;
                            if (proportion >= 0)
                                logger.info(proportion + "% completed");
                        }
                        lock.unlock();
                    }
                }
            };
        };

        this.logger.debug(this.statisticForTest.size() + " genes to be tested");
        for (String label: this.statisticForTest.keySet()) {
            Runnable runnable = task.runTask(label);
            threadPoolExecutor.execute(runnable);
        }
        try {
            countDown.await();
        } catch (InterruptedException ie) {
            this.logger.error("analysis interrupted");
            this.logger.error(ie.getMessage());
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
        this.logger.info("model ase test complete");
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
    
    public PValueParam getPValue(double predictAse, String geneId) {
        String distributionType = getDistributionType(geneId);
        double pvalue;
        if (distributionType == null) {
            distributionType = "unknown";
        }
        RecordParamDTO paramDTO = distributionType2GpdParam.getOrDefault(distributionType, null);
        GPDFunction gpdFunction = new GPDFunction(paramDTO.gpdGamma, paramDTO.gpdSigma);
        if (predictAse >= paramDTO.t) {
            pvalue = ((double) paramDTO.Nt / paramDTO.N) * (1 - gpdFunction.cdf(predictAse - paramDTO.t));
        } else {
            double count = 0;
            for (int i = 0; i < paramDTO.samplingMafs.size(); i++) {
                if (paramDTO.samplingMafs.get(i) >= predictAse) {
                    count++;
                }
            }
            pvalue = count / paramDTO.samplingMafs.size();
        }
        return new PValueParam(pvalue, distributionType);
    }
    
    private void bhRecalibrationOfEachGene() {
        this.logger.info("start recalibrating p values of hierarchical model");
        this.logger.info("sorting test result in order");
        // geneAsePvalue :   distributionType -> pvalue -> list(geneIdName, geneIdName)
        this.geneAseQValue = new ArrayList<>();
        for (String disType : geneAsePValue.keySet()) {
            HashMap<Double, ArrayList<String>> pvalue2GeneIdName = geneAsePValue.get(disType);
    
            ArrayList<Map.Entry<Double, ArrayList<String>>> sortedPVals = new ArrayList<>(pvalue2GeneIdName.entrySet());
            // sort p value from large to small
            sortedPVals.sort((o1, o2) -> o2.getKey().compareTo(o1.getKey()));
            
            int totalGenes = 0;
            for (Map.Entry<Double, ArrayList<String>> geneAse : pvalue2GeneIdName.entrySet()) {
                totalGenes += geneAse.getValue().size();
            }
            int rankage = totalGenes;
            double prevQValue = 1.0, qValue;
            String pValString, qValString;
    
            for (Map.Entry<Double, ArrayList<String>> entry: sortedPVals) {
                Double pVal = entry.getKey();
                ArrayList<String> samePValGenes = entry.getValue();
                // sort by SNV numbers
                HashMap<String, Integer> samePValGeneSNVs = new HashMap<>();
                for (String gene: samePValGenes) {
                    samePValGeneSNVs.put(gene, this.genesSNVs.get(gene));
                }
                // sort by major allele frequency(MAF)
                HashMap<String, Double> samePValGeneMajorAlleleFrequency = new HashMap<>();
                for (String gene: samePValGenes) {
                    samePValGeneMajorAlleleFrequency.put(gene, this.geneMajorAlleleFrequency.get(gene));
                }
        
                List<Map.Entry<String, Integer>> samePValGeneEntry = new ArrayList<>(samePValGeneSNVs.entrySet());
                samePValGeneEntry.sort((o1, o2) -> {
            
                    // first sort genes with same p value with their SNV number, then sort by MAF if get same SNV number
                    // both SNV number and MAF are sorted from large to small
                    String gene1 = o1.getKey(), gene2 = o2.getKey();
                    Double gene1MAF = samePValGeneMajorAlleleFrequency.get(gene1), gene2MAF = samePValGeneMajorAlleleFrequency.get(gene2);
                    Integer gene1SNVs = o1.getValue(), gene2SNVs = o2.getValue();
                    if (gene1SNVs.equals(gene2SNVs)) {
                        return gene2MAF.compareTo(gene1MAF);
                    } else
                        return gene2SNVs.compareTo(gene1SNVs);
                });
        
                for (Map.Entry<String, Integer> geneEntry: samePValGeneEntry) {
                    String geneIdName = geneEntry.getKey();
                    qValue = Math.min(prevQValue, pVal * totalGenes / rankage);
                    if (qValue - prevQValue < 0.00001)
                        prevQValue = qValue;
                    rankage--;
            
                    pValString = String.valueOf(pVal);
                    qValString = String.valueOf(qValue);
                    this.geneAseQValue.add(String.join("->", new String[]{geneIdName, pValString, qValString}));
                }
            }
    
        }
        this.logger.info("recalibration complete.");
    }

    /**
     * write the test result into file
     * statusCode : 0 -> asm algorithm call, need to specify the output file name  1 -> ase algorithm call
     */
    private void outputResult(int statusCode) {
        HashMap<String, String[]> finalRecords = new HashMap<>();
        BufferedWriter bfw = null;
        if (statusCode == 0) {
            String fileParentPath = new File(aseGeneFile).getParentFile().getAbsolutePath();
            aseGeneFile = fileParentPath + "/aseRes.txt";
        }
        try {
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(this.aseGeneFile)));
            String line, label, geneId, geneName, majorAlleleRecord, pVal, qVal;
            String[] info;
            int snvNum, majorCount, minorCount, majorBackground, minorBackground;
            double majorAlleleFrequency;
            bfw.write("geneId\tgeneName\tpValue\tqValue\tsnpNum\tmajor/minorAlleleReads\tmajorAlleleFrequency\tmajorAlleleBase\n");
            for (String record: this.geneAseQValue) {
                info = record.split("->");
                geneId = info[0];
                geneName = info[1];
                pVal = info[2];
                qVal = info[3];
                label = String.join("->", new String[]{geneId, geneName});
                majorCount = this.geneReadsCount.get(label).get("major");
                minorCount = this.geneReadsCount.get(label).get("minor");
                majorAlleleFrequency = this.geneMajorAlleleFrequency.get(label);
                majorAlleleRecord = this.geneMajorNucleotide.get(geneId);

                snvNum = this.genesSNVs.get(label);
                majorBackground = (this.geneBackgroundReads == null)? 0: this.geneBackgroundCount.get(label).get("major");
                minorBackground = (this.geneBackgroundReads == null)? 0: this.geneBackgroundCount.get(label).get("minor");

                finalRecords.put(label, new String[]{geneId, geneName, pVal, qVal, Integer.toString(snvNum),
                        Integer.toString(majorCount), Integer.toString(minorCount), Integer.toString(majorBackground),
                        Integer.toString(minorBackground), String.valueOf(majorAlleleFrequency), majorAlleleRecord});
            }

            List<Map.Entry<String, String[]>> records = new ArrayList<>(finalRecords.entrySet());
            // sort items with its q value
            records.sort((o1, o2) -> {
                String[] data1 = o1.getValue(), data2 = o2.getValue();
                Double q1 = Double.parseDouble(data1[3]), q2 = Double.parseDouble(data2[3]);
                if (!q1.equals(q2)) {
                    return q1.compareTo(q2);
                }
                // sort gene records with SNV numbers if has same q value
                Integer snvCount1 = Integer.parseInt(data1[4]), snvCount2 = Integer.parseInt(data2[4]);
                if (snvCount1 - snvCount2 != 0) {
                    return snvCount2.compareTo(snvCount1);
                }
                // sort gene records with MAF if has same q value
                Integer majorCount1 = Integer.parseInt(data1[5]), minorCount1 = Integer.parseInt(data1[6]),
                        majorCount2 = Integer.parseInt(data2[5]), minorCount2 = Integer.parseInt(data2[6]);
                Double maf1 = (double) majorCount1 / (double) (majorCount1 + minorCount1),
                        maf2 = (double) majorCount2 / (double) (majorCount2 + minorCount2);
                return maf2.compareTo(maf1);
            });

            for (Map.Entry<String, String[]> rec: records) {
                String[] data = rec.getValue();
                String[] lineInfo = new String[]{data[0], data[1], data[2], data[3], data[4], data[5] + "," + data[6], data[9], data[10]};
                line = String.join("\t", lineInfo);
                bfw.write(line);
                bfw.newLine();
            }
            this.logger.info("ase result file " + this.aseGeneFile);
        } catch (IOException ie) {
            this.logger.error(ie.getMessage());
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

    private Integer getSum(ArrayList<Integer> list) {
        Integer total = 0;
        for (Integer i: list)
            total += i;

        return total;
    }
    
    public Map<String, Double> getGeneMajorAlleleOddRatio() {
        return this.geneMajorAlleleOddRatio;
    }

    private String getString(LinkedList<String> list) {
        list.sort((o1, o2) -> {
            Integer pos1 = Integer.parseInt(o1.split(":")[0]);
            Integer pos2 = Integer.parseInt(o2.split(":")[0]);
            return pos2.compareTo(pos1);
        });

        String[] str = new String[list.size()];
        for (int i=0; i<list.size(); i++) {
            str[i] = list.get(i);
        }
    
        return String.join(";", str);
    }

    private static Logger initLog(String logHome) {
        System.setProperty("log_home", logHome);
        return Logger.getLogger(AseGeneDetection.class);
    }

    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("vcf", "vcf_file", true, "VCF format file generate by RNA-seq or MeRIP-seq data SNP calling process, required");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("g", "gtf", true, "GTF annotation file, required");
        option.setRequired(true);
        options.addOption(option);
        
        option = new Option("bam", "bam_files", true, "bamFiles, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);
    
        option = new Option("bai", "bai_files", true, "bamFile index file, required and separate with commas");
        option.setRequired(true);
        options.addOption(option);

        option = new Option("o", "output", true, "ASE gene test output directory. Optional, default .");
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

        option = new Option("h", "help", false, "help message of AseGeneDetection");
        option.setRequired(false);
        options.addOption(option);

        CommandLineParser parser = new DefaultParser();

        return parser.parse(options, args);
    }
}
