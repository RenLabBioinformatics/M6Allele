package HierarchicalBayesianAnalysis;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Arrays;
import java.util.HashMap;

public class HierarchicalBayesianModelAse {
    private final TauSampler ts;
    private double curTau, curGlobalLOR;
    private final int samplingTime;
    private final int burnIn;
    private boolean allZero = true;
    // vector of y, sigma, p(y_j|theta_j), theta
    private double[] observeLogOddRatio;
    private double[] variances;
    private final double[] singleASELORMean;
    // u sampling list
    private double[] samplingGlobalLORs = null;
    private final int[] majorAlleleReads;
    private final int[] minorAlleleReads;
    private final int[] majorAlleleBackground;
    private final int[] minorAlleleBackground;
    private final double df;
    private final double scaleParam;
    
    public HierarchicalBayesianModelAse(double df, double scaleParam, int samplingTime, int burnIn,
                                        int[] majorAlleleReads, int[] minorAlleleReads,
                                        int[] majorAlleleBackground, int[] minorAlleleBackground) {
        this.ts = new TauSampler(df, scaleParam);
        this.samplingTime = samplingTime;
        this.burnIn = burnIn;
        this.majorAlleleReads = majorAlleleReads;
        this.minorAlleleReads = minorAlleleReads;
        this.singleASELORMean = new double[minorAlleleReads.length];
        this.majorAlleleBackground = majorAlleleBackground;
        this.minorAlleleBackground = minorAlleleBackground;
        this.df = df;
        this.scaleParam = scaleParam;
    }

    /**
     * test for ASE or ASM significance. Return significant p value
     */
    public void testSignificant() {
        this.checkLORVectors();
        if (this.allZero)
            return;
        this.initializer();
        samplingGlobalLORs = KeanuSampling.sampleUseMh(df, scaleParam, variances, observeLogOddRatio,
                curGlobalLOR, curTau, singleASELORMean, samplingTime, burnIn);
    }
    
    private void checkLORVectors() {
        int idx = 0;
        while (idx < this.majorAlleleReads.length && this.allZero) {
            if (this.majorAlleleReads[idx] != this.minorAlleleReads[idx]) {
                this.allZero = false;
                return;
            }
            idx++;
        }
    }

    /**
     * hierarchical model parameters initialization
     */
    private void initializer() {
        // init LOR(vector y) and variance(vector sigma)
        HashMap<String, double[]> initLORAndVar = this.getObserveData();
        this.observeLogOddRatio = initLORAndVar.get("LOR");
        this.variances = initLORAndVar.get("VAR");

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
    }

    /**
     * calculate LOR and variance via observed data
     * @return {"LOR": [log odd ratios], "VAR": [variances]}
     */
    private HashMap<String, double[]> getObserveData() {
        OddRatioCalcAse orc = new OddRatioCalcAse(this.majorAlleleReads, this.minorAlleleReads,
                                            this.majorAlleleBackground, this.minorAlleleBackground);
        return orc.getLogOddRatio();
    }

    /**
     * quantify the global log odd ratio with sampling median
     * @return gene logarithm odd ratio
     */
    public double quantifyGeneLOR() {
        double median;
        if (this.allZero)
            median = 0;
        else {
            double[] sortedSamplingValue = Arrays.stream(this.samplingGlobalLORs).sorted().toArray();
            int medianIdx = sortedSamplingValue.length / 2;
            if (sortedSamplingValue.length % 2 == 0)
                median = (sortedSamplingValue[medianIdx] + sortedSamplingValue[ medianIdx - 1]) / 2;
            else
                median = sortedSamplingValue[medianIdx];

            this.samplingGlobalLORs = null;
        }

        return median;
    }
    
}
