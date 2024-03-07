package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class HierarchicalBayesianModelAsm {
    private final TauSampler ts;
    private double curTau;
    private double curGlobalLOR;
    private final int samplingTime;
    private final int burnIn;
    // vector of y, sigma, p(y_j|theta_j), theta
    private double[] observeLogOddRatio;
    private double[] variances;
    private final double[] singleASELORMean;
    private List<String> validSnpInfo;
    // u sampling list
    private double[] samplingGlobalLORs = null;
    // [majorNc:cnt:bam1:snpPos, minorNc:cnt:bam1:snpPos, majorNc:cnt:bam2:snpPos,...]
    private final String[] majorAlleleReads;
    private final String[] minorAlleleReads;
    private final double GeneExpressOddRatio;
    private final HashMap<Integer, String[]> aseSnpPos2BaseCount;
    private final double df;
    private final double scaleParam;
    private Integer snpNum;
    /**
     * Constructor
     */
    public HierarchicalBayesianModelAsm(double df, double scaleParam, int samplingTime, int burnIn,
                                        String[] majorAlleleReads, String[] minorAlleleReads, double GeneExpressOddRatio,
                                        HashMap<Integer, String[]> aseSnpPos2BaseCount) {
        this.ts = new TauSampler(df, scaleParam);
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.majorAlleleReads = majorAlleleReads;
        this.minorAlleleReads = minorAlleleReads;
        this.singleASELORMean = new double[minorAlleleReads.length];
        this.GeneExpressOddRatio = GeneExpressOddRatio;
        this.aseSnpPos2BaseCount = aseSnpPos2BaseCount;
        this.df = df;
        this.scaleParam = scaleParam;
    }
    
    public Integer getSnpNum() {
        return snpNum;
    }
    
    public List<String> getValidSnpInfo() {
        return validSnpInfo;
    }
    
    /**
     * mcmc sampling
     */
    public boolean testSignificant() {
        boolean flag = this.initializer();
        if (!flag) {
            return false;
        }
        this.snpNum = this.observeLogOddRatio.length;
        samplingGlobalLORs = KeanuSampling.sampleUseMh(df, scaleParam, variances, observeLogOddRatio,
                curGlobalLOR, curTau, singleASELORMean, samplingTime, burnIn);
        return true;
    }

    /**
     * hierarchical model parameters initialization
     */
    private boolean initializer() {
        // init LOR(vector y) and variance(vector sigma)
        HashMap<String, List<String>> initLORAndVar = this.getObserveData();
        this.observeLogOddRatio = initLORAndVar.get("LOR").stream().mapToDouble(Double::parseDouble).toArray();
        this.variances = initLORAndVar.get("VAR").stream().mapToDouble(Double::parseDouble).toArray();
        if (this.observeLogOddRatio.length == 0 || this.variances.length == 0) {
            return false;
        }
        validSnpInfo = initLORAndVar.get("majorInfo");

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

    /**
     * calculate LOR and variance via observed data
     * @return {"LOR": [log odd ratios], "VAR": [variances]}
     */
    private HashMap<String, List<String>> getObserveData() {
        OddRatioCalcAsm orc = new OddRatioCalcAsm(this.majorAlleleReads, this.minorAlleleReads, this.GeneExpressOddRatio,
                this.aseSnpPos2BaseCount);
        return orc.getLogOddRatio();
    }
    
    /**
     * quantify the global log odd ratio with sampling median
     * @return gene logarithm odd ratio
     */
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
