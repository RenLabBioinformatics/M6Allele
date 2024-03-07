package DifferentiationAnalysis;

import HierarchicalBayesianAnalysis.KeanuSampling;
import HierarchicalBayesianAnalysis.TauSampler;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;

/**
 * @author: tanglin
 * @date: 2023/09/22 15:35
 * @Description:
 */
public class SpecificAsmBayesianModel {
    
    private final TauSampler ts;
    private double curTau;
    private double curGlobalLOR;
    private final int samplingTime;
    private final int burnIn;
    // vector of y, sigma, p(y_j|theta_j), theta
    private double[] observeLogOddRatio, variances, singleASELORMean;
    // u sampling list
    private double[] samplingGlobalLORs = null;
    // [majorNc:cnt:bam1:snpPos, minorNc:cnt:bam1:snpPos, majorNc:cnt:bam2:snpPos,...]
    private final String[] sample1MajorAlleleReads;
    private final String[] sample1MinorAlleleReads;
    private final String[] sample2MajorAlleleReads;
    private final String[] sample2MinorAlleleReads;
    private final double sample1GeneExpressOddRatio;
    private final double sample2GeneExpressOddRatio;
    // snp -> [majorNc:cnt:bam1, minorNc:cnt:bam1, ...]
    private final Map<Integer, String[]> sample1AseSnpPos2BaseCount;
    private final Map<Integer, String[]> sample2AseSnpPos2BaseCount;
    private final double df;
    private final double scaleParam;
    private int snpNum;
    private List<String> validSnpInfo = new ArrayList<>();
    
    public SpecificAsmBayesianModel(int samplingTime, int burnIn,
                                    String[] sample1MajorAlleleReads, String[] sample1MinorAlleleReads,
                                    String[] sample2MajorAlleleReads, String[] sample2MinorAlleleReads,
                                    double sample1GeneExpressOddRatio, double sample2GeneExpressOddRatio,
                                    Map<Integer, String[]> sample1AseSnpPos2BaseCount,
                                    Map<Integer, String[]> sample2AseSnpPos2BaseCount,
                                    double df, double scaleParam) {
        this.ts = new TauSampler(df, scaleParam);
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.sample1MajorAlleleReads = sample1MajorAlleleReads;
        this.sample1MinorAlleleReads = sample1MinorAlleleReads;
        this.sample2MajorAlleleReads = sample2MajorAlleleReads;
        this.sample2MinorAlleleReads = sample2MinorAlleleReads;
        this.sample1GeneExpressOddRatio = sample1GeneExpressOddRatio;
        this.sample2GeneExpressOddRatio = sample2GeneExpressOddRatio;
        this.sample1AseSnpPos2BaseCount = sample1AseSnpPos2BaseCount;
        this.sample2AseSnpPos2BaseCount = sample2AseSnpPos2BaseCount;
        this.df = df;
        this.scaleParam = scaleParam;
    }
    
    public int getSnpNum() {
        return snpNum;
    }
    
    public List<String> getValidSnpInfo() {
        return validSnpInfo;
    }
    
    public double testSignificant() {
        boolean flag = initializer();
        if (!flag) {
            return Double.NaN;
        }
        snpNum = observeLogOddRatio.length;
        samplingGlobalLORs = KeanuSampling.sampleUseMh(df, scaleParam, variances, observeLogOddRatio,
                curGlobalLOR, curTau, singleASELORMean, samplingTime, burnIn);
        double positiveLOR = 0, negativeLOR = 0;
        for (double globalLOR: this.samplingGlobalLORs) {
            if (globalLOR <= 0.0000001)
                negativeLOR++;
            else
                positiveLOR++;
        }
        return 2 * Math.min(positiveLOR, negativeLOR) / this.samplingTime;
    }
    
    
    private boolean initializer() {
        Map<String, List<String>> initLorVar = this.getObserveData();
        this.observeLogOddRatio = initLorVar.get("LOR").stream().mapToDouble(Double::parseDouble).toArray();
        this.variances = initLorVar.get("VAR").stream().mapToDouble(Double::parseDouble).toArray();
        this.singleASELORMean = new double[observeLogOddRatio.length];
        if (this.observeLogOddRatio.length == 0 || this.variances.length == 0) {
            return false;
        }
        validSnpInfo = initLorVar.get("validSnpInfo");
        // randomly init value tau from priority distribution
        this.curTau = this.ts.randomInit();
    
        // init u with known tau, y and sigma
        double std, miuDenominator = 0, miuNumerator = 0;
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            std = this.variances[i] + Math.pow(this.curTau, 2);
            miuDenominator += 1.0 / std;
            miuNumerator += 1.0 / std * this.observeLogOddRatio[i];
        }
        double globalLORMean = miuNumerator / miuDenominator; // u^
        double globalLORSigma = 1.0 / miuDenominator;         // Vu
        NormalDistribution nd = new NormalDistribution(globalLORMean, Math.sqrt(globalLORSigma));
        this.curGlobalLOR = nd.sample();
    
        // init vector theta with known u, tau, y and sigma
        for (int i = 0; i < this.observeLogOddRatio.length; i++) {
            double mean = (1.0 / this.variances[i] * this.observeLogOddRatio[i] + 1.0 / Math.pow(this.curTau, 2) * this.curGlobalLOR) / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2)); // theta^_j
            double var = 1.0 / (1.0 / this.variances[i] + 1.0 / Math.pow(this.curTau, 2));  // V_j
            nd = new NormalDistribution(mean, Math.sqrt(var));
            double lor = nd.sample();
            this.singleASELORMean[i] = lor;
        }
        
        return true;
    }
    
    private Map<String, List<String>> getObserveData() {
        List<String> logOR = new ArrayList<>();
        List<String> variance = new ArrayList<>();
        List<String> validSnpInfo = new ArrayList<>();
        String sample1MajorReads, sample1MinorReads;
        String sample2MajorReads, sample2MinorReads;
        double logOddRatio, var;
        for (int i = 0; i < this.sample1MajorAlleleReads.length; i++) {
            sample1MajorReads = sample1MajorAlleleReads[i];
            sample1MinorReads = sample1MinorAlleleReads[i];
            sample2MajorReads = sample2MajorAlleleReads[i];
            sample2MinorReads = sample2MinorAlleleReads[i];
            
            logOddRatio = this.calculateLogOddRatio(sample1MajorReads, sample1MinorReads,
                    sample2MajorReads, sample2MinorReads);
            if (Double.isNaN(logOddRatio)) {
                continue;
            }
            var = calculateOneSampleVariance(sample1MajorReads, sample1MinorReads,
                    sample1AseSnpPos2BaseCount, sample1GeneExpressOddRatio) +
                    calculateOneSampleVariance(sample2MajorReads, sample2MinorReads,
                            sample2AseSnpPos2BaseCount, sample2GeneExpressOddRatio);
            logOR.add(String.valueOf(logOddRatio));
            variance.add(String.valueOf(var));
            validSnpInfo.add(sample1MajorReads);
            validSnpInfo.add(sample1MinorReads);
        }
        HashMap<String, List<String>> calcRes = new HashMap<>();
        calcRes.put("LOR", logOR);
        calcRes.put("VAR", variance);
        calcRes.put("validSnpInfo", validSnpInfo);
        
        return calcRes;
    }
    
    private double calculateLogOddRatio(String sample1MajorAlleleRead, String sample1MinorAlleleRead,
                                        String sample2MajorAlleleRead, String sample2MinorAlleleRead) {
    
        double sample1OddRatio = calculateOneSampleLogOddRatio(sample1MajorAlleleRead,
                sample1MinorAlleleRead, sample1AseSnpPos2BaseCount, sample2GeneExpressOddRatio);
        if (Double.isNaN(sample1OddRatio)) {
            return Double.NaN;
        }
        double sample2OddRatio = calculateOneSampleLogOddRatio(sample2MajorAlleleRead,
                sample2MinorAlleleRead, sample2AseSnpPos2BaseCount, sample2GeneExpressOddRatio);
        if (Double.isNaN(sample2OddRatio)) {
            return Double.NaN;
        }
        return sample1OddRatio - sample2OddRatio;
    
    }
    
    private double calculateOneSampleLogOddRatio(String majorAlleleReads, String minorAlleleReads,
                                                 Map<Integer, String[]> aseSnpPos2BaseCount,
                                                 double geneExpressOddRatio) {
        double oddRatio = -1;
        try {
            // majorNc:cnt:bam:snpPos
            String[] majorInfos = majorAlleleReads.split(":");
            String[] minorInfos = minorAlleleReads.split(":");
            int snpPos = Integer.parseInt(majorInfos[3]);
            String[] allele2Count = aseSnpPos2BaseCount.get(snpPos);
            if (allele2Count == null || allele2Count.length == 0) {
                return Double.NaN;
            }
            // correspondAseSnpInfo:  majorNc:count:bam1, minorNc:count:bam1
            String[] correspondAseSnpInfo = getCorrespondAseSnpInfo(allele2Count, majorInfos[2]);
            if (correspondAseSnpInfo == null) {
                return Double.NaN;
            }
            String aseMajor = correspondAseSnpInfo[0].split(":")[0];
            String aseMinor = correspondAseSnpInfo[1].split(":")[0];
            String asmMajor = majorInfos[0];
            String asmMinor = minorInfos[0];
            if (aseMajor.equalsIgnoreCase(asmMajor) && aseMinor.equalsIgnoreCase(asmMinor)) {
                oddRatio = (Double.parseDouble(majorInfos[1]) + 0.1) / (Double.parseDouble(minorInfos[1]) + 0.1) / geneExpressOddRatio;
            } else if (aseMajor.equalsIgnoreCase(asmMinor) && aseMinor.equalsIgnoreCase(asmMajor)) {
                oddRatio = (Double.parseDouble(minorInfos[1]) + 0.1) / (Double.parseDouble(majorInfos[1]) + 0.1) / geneExpressOddRatio;
            }
        } catch (Exception e) {
            return Double.NaN;
        }
        return Math.log(oddRatio);
    }
    
    private double calculateOneSampleVariance(String majorAlleleReads, String minorAlleleReads,
                                              Map<Integer, String[]> aseSnpPos2BaseCount,
                                              double geneExpressOddRatio) {
        // majorNc:cnt:bam:snpPos
        String[] majorInfos = majorAlleleReads.split(":");
        String[] minorInfos = minorAlleleReads.split(":");
        int asmMajorCnt = Integer.parseInt(majorInfos[1]);
        int asmMinorCnt = Integer.parseInt(minorInfos[1]);
        String bamFileName = majorInfos[2];
        int snpPos = Integer.parseInt(majorInfos[3]);
        String[] allele2Count = aseSnpPos2BaseCount.get(snpPos);
        // majorNc:cnt:bam, minorNc:cnt:bam
        String[] correspondAseSnpInfo = getCorrespondAseSnpInfo(allele2Count, bamFileName);
        assert correspondAseSnpInfo != null;
        String aseMajorNc = correspondAseSnpInfo[0].split(":")[0];
        String aseMinorNc = correspondAseSnpInfo[1].split(":")[0];
        if (majorInfos[0].equals(aseMajorNc) && minorInfos[0].equals(aseMinorNc)) {
            return 1.0 / ((double) asmMajorCnt / asmMinorCnt) + 1 / geneExpressOddRatio;
        } else {
            return 1.0 / ((double) asmMinorCnt / asmMajorCnt) + 1 / geneExpressOddRatio;
        }
    }
    /**
     *
     * @param aseSnpInfo  [majorNc:majorCount:bam1, minorNc:minorCount:bam1, ...]
     * @param bamFileName bam1 || bam2...
     * @return
     */
    private String[] getCorrespondAseSnpInfo(String[] aseSnpInfo, String bamFileName) {
        for (int i = 0; i < aseSnpInfo.length; i++) {
            String snpInfo = aseSnpInfo[i];
            if (snpInfo.contains(bamFileName)) {
                return new String[] {aseSnpInfo[i], aseSnpInfo[i + 1]};
            }
        }
        return null;
    }
    
    public double quantifyGeneLOR() {
        double median;
        double[] sortedSamplingValue = Arrays.stream(this.samplingGlobalLORs).sorted().toArray();
        int medianIdx = sortedSamplingValue.length / 2;
        if (sortedSamplingValue.length % 2 == 0) {
            median = (sortedSamplingValue[medianIdx] + sortedSamplingValue[medianIdx + 1]) / 2;
        } else {
            median = sortedSamplingValue[medianIdx];
        }
        this.samplingGlobalLORs = null;
        return median;
    }
    
    
}
