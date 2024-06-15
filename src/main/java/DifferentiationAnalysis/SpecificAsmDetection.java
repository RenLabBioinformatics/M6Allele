package DifferentiationAnalysis;

import HierarchicalBayesianAnalysis.AseGeneDetection;
import HierarchicalBayesianAnalysis.AsmPeakDetection;
import HierarchicalBayesianAnalysis.PValueParam;
import HierarchicalBayesianAnalysis.RunTest;
import Threshold.GPD.GPDFunction;
import Threshold.RecordParamDTO;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;

/**
 * @author: tanglin
 * @date: 2023/09/22 15:35
 * @Description:
 */
public class SpecificAsmDetection {
    
    // chr:peakStart:peakEnd:geneId ->
    // {pos1: [majorNc:cnt:bam1, minorNc:cnt:bam1, majorNc:cnt:bam2, minorNc:cnt:bam2...], pos2:...}
    private final HashMap<String, HashMap<String, List<String>>> sample1PeakSnpReadsCount;
    private final HashMap<String, HashMap<String, List<String>>> sample2PeakSnpReadsCount;
    // geneId->geneName oddRadio
    private final Map<String, Double> sample1GeneOddRadio;
    private final Map<String, Double> sample2GeneOddRadio;
    // {"geneId->geneName": {pos1: [majorAllele:count:bam1, minorAllele:count:bam1],...}, ...}
    private final HashMap<String, HashMap<Integer, String[]>> sample1GeneAlleleReads;
    private final HashMap<String, HashMap<Integer, String[]>> sample2GeneAlleleReads;
    // [chr:peakStart:peakEnd:geneId, ...]
    private final Set<String> needResamplingPeakLabels;
    private Map<String, List<String[]>> statisticForTest = new HashMap<>();
    private final Logger log;
    private final int threadNumber;
    private final HashMap<String, HashMap<String, String>> chr2GeneId2Name;
    private final int sampleTimes;
    private final int burnIn;
    private final Map<Double, List<String>> sameSpecificPValue = new HashMap<>();
    //private final Map<String, Double> peakMajorLogOddRatio = new HashMap<>();
    private final Map<String, Double> peakMajorMAF = new HashMap<>();
    private ReentrantLock lock;
    private final Map<String, Integer> peakLabel2SnpSize = new HashMap<>();
    private final List<String> pairedCompareQValue = new ArrayList<>();
    private final List<String> discardPeak = new CopyOnWriteArrayList<>();
    private final Map<String, List<String>> peakLabel2ValidSnpInfo = new HashMap<>();
    private final Map<String, LinkedHashSet<String>> peakMajorNc = new HashMap<>();
    public Map<String, RecordParamDTO> distributionType2GpdParam = new HashMap<>();
    public Map<String, List<String>> distributionType2GeneIds = new HashMap<>();
    
    public SpecificAsmDetection(HashMap<String, HashMap<String, List<String>>> sample1PeakSnpReadsCount,
                                HashMap<String, HashMap<String, List<String>>> sample2PeakSnpReadsCount,
                                Map<String, Double> sample1GeneOddRadio, Map<String, Double> sample2GeneOddRadio,
                                HashMap<String, HashMap<Integer, String[]>> sample1GeneAlleleReads,
                                HashMap<String, HashMap<Integer, String[]>> sample2GeneAlleleReads,
                                Set<String> needResamplingPeakLabels, Logger logger, int threadNumber,
                                HashMap<String, HashMap<String, String>> geneNames, int sampleTimes, int burnIn) {
        this.sample1PeakSnpReadsCount = sample1PeakSnpReadsCount;
        this.sample2PeakSnpReadsCount = sample2PeakSnpReadsCount;
        this.sample1GeneOddRadio = sample1GeneOddRadio;
        this.sample2GeneOddRadio = sample2GeneOddRadio;
        this.sample1GeneAlleleReads = sample1GeneAlleleReads;
        this.sample2GeneAlleleReads = sample2GeneAlleleReads;
        this.needResamplingPeakLabels = needResamplingPeakLabels;
        this.log = logger;
        this.threadNumber = threadNumber;
        this.chr2GeneId2Name = geneNames;
        this.sampleTimes = sampleTimes;
        this.burnIn = burnIn;
        getGpdFunction();
        getDistributionType2GeneIds();
    }
    
    public void asmPairedCompare() {
        this.dataPreparation();
        this.hierarchicalModelTest();
        this.recordValidSnpInfo();
        this.bhRecalibrationOfEachPeak();
    }
    
    public List<String> getPairedCompareQValue() {
        return pairedCompareQValue;
    }
    
    public List<String> getDiscardPeak() {
        return discardPeak;
    }
    
    public Map<String, LinkedHashSet<String>> getPeakMajorNc() {
        return peakMajorNc;
    }
    
    public void hierarchicalModelTest() {
        this.log.info("start sampling paired sample");
        ExecutorService threadPoolExecutor = Executors.newFixedThreadPool(this.threadNumber);
        lock = new ReentrantLock();
        CountDownLatch countDown = new CountDownLatch(statisticForTest.size());
        long tenPercent = Math.round(statisticForTest.size() * 0.1);
        RunTest task = (String label) -> {
            return new Runnable() {
                @Override
                public void run() {
                    List<String[]> statistic;
                    double df = 5, scaleParam = 10, s1GeneExpressOddRation, s2GeneExpressOddRatio;
                    String[] sample1MajorPeakReads, sample1MinorPeakReads,
                            sample2MajorPeakReads, sample2MinorPeakReads;
                    SpecificAsmBayesianModel model;
                    try {
                        statistic = statisticForTest.get(label);
                        sample1MajorPeakReads = statistic.get(0);
                        sample1MinorPeakReads = statistic.get(1);
                        sample2MajorPeakReads = statistic.get(2);
                        sample2MinorPeakReads = statistic.get(3);
                        // chr:peakStart:peakEnd:geneId
                        String[] labelArr = label.split(":");
                        String geneId = labelArr[3];
                        String geneName = chr2GeneId2Name.get(labelArr[0]).getOrDefault(geneId, "unknown");
                        String name = String.join("->", new String[]{geneId, geneName});
                        
                        if (sample1GeneOddRadio.get(name) == null || sample2GeneOddRadio.get(name) == null) {
                            return;
                        } else {
                            s1GeneExpressOddRation = sample1GeneOddRadio.get(name);
                            s2GeneExpressOddRatio = sample2GeneOddRadio.get(name);
                        }
                        // snpPos -> [majorAllele:Count:bamFile1, minorAllele:count:bamFile1,
                        // majorAllele:count:bam2, minorAllele:count:bam2....] 多个bam文件
                        HashMap<Integer, String[]> sample1SnpPos2BaseCount = sample1GeneAlleleReads.get(name);
                        HashMap<Integer, String[]> sample2SnpPos2BaseCount = sample2GeneAlleleReads.get(name);
                        // get p value via hierarchical model
                        double pVal, peakOddRatio, peakMAF;
                        model = new SpecificAsmBayesianModel(sampleTimes, burnIn, sample1MajorPeakReads, sample1MinorPeakReads,
                                sample2MajorPeakReads, sample2MinorPeakReads, s1GeneExpressOddRation, s2GeneExpressOddRatio,
                                sample1SnpPos2BaseCount, sample2SnpPos2BaseCount, df, scaleParam);
                        pVal = model.testSignificant();
                        
                        if (Double.isNaN(pVal)) {
                            log.debug(label + "peak is discard");
                            discardPeak.add(label);
                            return;
                        }
                        peakOddRatio = Math.exp(model.quantifyGeneLOR());
                        peakMAF = Math.min(1.0, peakOddRatio / (peakOddRatio + 1));
                        PValueParam pValue = getPValue(peakMAF, geneId);
                        pVal = pValue.pvalue;
                        lock.lock();
                        
                        peakLabel2SnpSize.put(label, model.getSnpNum());
                        peakLabel2ValidSnpInfo.put(label, model.getValidSnpInfo());
                        List<String> samePValPeaks = sameSpecificPValue.getOrDefault(pVal, new ArrayList<>());
                        samePValPeaks.add(label);
                        sameSpecificPValue.put(pVal, samePValPeaks);
                        peakMajorMAF.put(label, peakMAF);
                        
                    } catch (Exception e) {
                        log.error("sampling paired sample error occurs on record " + label);
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
        
        log.debug(statisticForTest.size() + " m6a peaks to be paired compare");
        for (String label : statisticForTest.keySet()) {
            Runnable runnable = task.runTask(label);
            threadPoolExecutor.execute(runnable);
        }
        try {
            countDown.await();
        } catch (InterruptedException e) {
            this.log.error("paired compare analysis success");
            this.log.error(e.getMessage());
        } finally {
            this.statisticForTest = null;
            try {
                threadPoolExecutor.shutdown();
                if (!threadPoolExecutor.awaitTermination(1000, TimeUnit.MINUTES)) {
                    threadPoolExecutor.shutdown();
                }
            } catch (InterruptedException e) {
                threadPoolExecutor.shutdown();
            }
            this.lock = null;
        }
    }
    
    public void dataPreparation() {
        for (String label : needResamplingPeakLabels) {
            HashMap<String, List<String>> sample1PeakSnpReads = sample1PeakSnpReadsCount.get(label);
            List<String> sample1MajorPeakReads = new ArrayList<>();
            List<String> sample1MinorPeakReads = new ArrayList<>();
            processPeakReads(sample1MajorPeakReads, sample1MinorPeakReads, sample1PeakSnpReads);
            HashMap<String, List<String>> sample2PeakSnpReads = sample2PeakSnpReadsCount.get(label);
            List<String> sample2MajorPeakReads = new ArrayList<>();
            List<String> sample2MinorPeakReads = new ArrayList<>();
            processPeakReads(sample2MajorPeakReads, sample2MinorPeakReads, sample2PeakSnpReads);
            String[] sample1MajorPeakReadsArr = sample1MajorPeakReads.toArray(new String[0]);
            String[] sample1MinorPeakReadsArr = sample1MinorPeakReads.toArray(new String[0]);
            String[] sample2MajorPeakReadsArr = sample2MajorPeakReads.toArray(new String[0]);
            String[] sample2MinorPeakReadsArr = sample2MinorPeakReads.toArray(new String[0]);
            List<String[]> records = new ArrayList<>();
            records.add(sample1MajorPeakReadsArr);
            records.add(sample1MinorPeakReadsArr);
            records.add(sample2MajorPeakReadsArr);
            records.add(sample2MinorPeakReadsArr);
            this.statisticForTest.put(label, records);
        }
    }
    private void processPeakReads(List<String> majorPeakReads, List<String> minorPeakReads,
                                  HashMap<String, List<String>> peakSnpReads) {
        for (String position : peakSnpReads.keySet()) {
            // [majorNc:cnt:bam1, minorNc:cnt:bam1, ...]
            List<String> snpInfos = peakSnpReads.get(position);
            for (int i = 0; i < snpInfos.size(); i += 2) {
                majorPeakReads.add(snpInfos.get(i) + ":" + position);
                minorPeakReads.add(snpInfos.get(i + 1) + ":" + position);
            }
        }
    }
    
    private void recordValidSnpInfo() {
        // peakLabel2ValidSnpInfo  chr:peakStart:peakEnd:geneId -> [majorNc:cnt:bam1:snpPos, minorNc:cnt:bam1:snpPos, ...]
        for (String label : peakLabel2ValidSnpInfo.keySet()) {
            List<String> snpInfos = peakLabel2ValidSnpInfo.get(label);
            LinkedHashSet<String> snpPosMajorNc = new LinkedHashSet<>();
            for (int i = 0; i < snpInfos.size(); i += 2) {
                String[] majorInfos = snpInfos.get(i).split(":");
                snpPosMajorNc.add(String.join(":", new String[]{majorInfos[3], majorInfos[0]}));
            }
            peakMajorNc.put(label, snpPosMajorNc);
        }
    }
    
    /**
     * get gpd function to compute threshold
     */
    public void getGpdFunction() {
        InputStream paramInputStream = AseGeneDetection.class.getClassLoader().getResourceAsStream("gpdParameterTwoSample.txt");
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
    
    private void bhRecalibrationOfEachPeak() {
        ArrayList<Map.Entry<Double, List<String>>> sortedPVals = new ArrayList<>(sameSpecificPValue.entrySet());
        sortedPVals.sort((a, b) -> b.getKey().compareTo(a.getKey()));
        int totalPeak = 0;
        for (Map.Entry<Double, List<String>> sortedPVal : sortedPVals) {
            totalPeak += sortedPVal.getValue().size();
        }
        int rankage = totalPeak;
        double preQValue = 1, qValue;
        String pValString, qValString;
        for (Map.Entry<Double, List<String>> entry : sortedPVals) {
            double pVal = entry.getKey();
            List<String> samePValPeaks = entry.getValue();
            Map<String, Integer> samePValPeaksSNPs = new HashMap<>();
            for (String peak : samePValPeaks) {
                samePValPeaksSNPs.put(peak, peakLabel2SnpSize.get(peak));
            }
            HashMap<String, Double> samePValPeakMajorAlleleFrequency = new HashMap<>();
            for (String peakLabel : samePValPeaks) {
                samePValPeakMajorAlleleFrequency.put(peakLabel, this.peakMajorMAF.get(peakLabel));
            }
            ArrayList<Map.Entry<String, Integer>> samePValPeakEntry = new ArrayList<>(samePValPeaksSNPs.entrySet());
            samePValPeakEntry.sort((o1, o2) -> {
                Integer peak1SNPs = o1.getValue();
                Integer peak2SNPs = o2.getValue();
                Double peak1MAF = samePValPeakMajorAlleleFrequency.get(o1.getKey());
                Double peak2MAF = samePValPeakMajorAlleleFrequency.get(o2.getKey());
                if (peak1SNPs.equals(peak2SNPs)) {
                    return peak2MAF.compareTo(peak1MAF);
                } else {
                    return peak2SNPs.compareTo(peak1SNPs);
                }
            });
            for (Map.Entry<String, Integer> peakEntry : samePValPeakEntry) {
                String peak = peakEntry.getKey();
                qValue = Math.min(preQValue, pVal * totalPeak / rankage);
                if (qValue - preQValue < 0.00001) {
                    preQValue = qValue;
                }
                rankage--;
                pValString = String.valueOf(pVal);
                qValString = String.valueOf(qValue);
                this.pairedCompareQValue.add(String.join("->",
                        new String[]{peak, pValString, qValString}));
            }
        }
    }
    
}
