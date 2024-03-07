package HierarchicalBayesianAnalysis;

import java.util.HashMap;

/**
 * calculate the LOR and its variance of SNV sites locate in the range of a gene or m6A signal
 */
public class OddRatioCalcAse {
    private final int[] majorSNPReads;
    private final int[] minorSNPReads;
    private final int[] majorAlleleBackground;
    private final int[] minorAlleleBackground;

    /**
     * Constructor
     * @param majorSNPReads MeRIP-seq INPUT data major allele reads count
     * @param minorSNPReads MeRIP-seq INPUT data minor allele reads count
     * @param majorAlleleBackground background major allele reads count
     * @param minorAlleleBackground background minor allele reads count
     */
    public OddRatioCalcAse(int[] majorSNPReads, int[] minorSNPReads, int[] majorAlleleBackground, int[] minorAlleleBackground) {
        this.majorSNPReads = majorSNPReads;
        this.minorSNPReads = minorSNPReads;
        this.majorAlleleBackground = majorAlleleBackground;
        this.minorAlleleBackground = minorAlleleBackground;
    }

    /**
     * calculate major allele LOR and the variance
     * @return {"LOR": [log orr ratios], "VAR": [variances]}
     */
    public HashMap<String, double[]> getLogOddRatio() {
        double[] logOR = new double[this.majorSNPReads.length];
        double[] variance = new double[this.majorSNPReads.length];
        int majorReads, minorReads, majorBackground, minorBackground;
        double logOddratio, var;
        for (int i = 0; i < this.majorSNPReads.length; i++) {
            majorReads = this.majorSNPReads[i];
            minorReads = this.minorSNPReads[i];
            majorBackground = this.majorAlleleBackground[i];
            minorBackground = this.minorAlleleBackground[i];
            
            logOddratio = this.calculateLogOddRatio(majorReads, minorReads, majorBackground, minorBackground);
            var = this.calculateVariance(majorReads, minorReads, majorBackground, minorBackground);
            logOR[i] = logOddratio;
            variance[i] = var;
        }
        HashMap<String, double[]> calcRes = new HashMap<>();
        calcRes.put("LOR", logOR);
        calcRes.put("VAR", variance);

        return calcRes;
    }

    /**
     * formula of LOR calculation
     *      yi = ln[y_major/(total-y_major)] - ln[0.5*total/(total-0.5*total)]
     *         = ln[y_major/(total-y_major)]
     * @param majorAlleleReads MeRIP-seq INPUT data major allele reads count
     * @param minorAlleleReads MeRIP-seq INPUT data minor allele reads count
     * @param majorBackground background major allele reads count
     * @param minorBackground background minor allele reads count
     * @return LOR of a SNV site
     */
    private double calculateLogOddRatio(int majorAlleleReads, int minorAlleleReads, int majorBackground, int minorBackground) {
        double oddRatio;
        oddRatio = (majorAlleleReads + 0.1) / (minorAlleleReads + 0.1) / ((majorBackground + 0.1) / (minorBackground + 0.1));
        return Math.log(oddRatio);
    }

    /**
     * formula of LOR variance calculation
     *      s = 1/y_major + 1/(total-y_major) + 1/majorBackground + 1/minorBackground
     * @param majorAlleleReads MeRIP-seq INPUT data major allele reads count
     * @param minorAlleleReads MeRIP-seq INPUT data minor allele reads count
     * @param majorBackground background major allele reads count
     * @param minorBackground background minor allele reads count
     * @return LOR variance of a SNV site
     */
    private double calculateVariance(int majorAlleleReads, int minorAlleleReads, int majorBackground, int minorBackground) {
        return  1.0 / majorAlleleReads + 1.0 / minorAlleleReads + 1.0 / majorBackground + 1.0 / minorBackground;
    }
}
