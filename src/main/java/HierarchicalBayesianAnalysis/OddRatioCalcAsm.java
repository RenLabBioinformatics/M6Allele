package HierarchicalBayesianAnalysis;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class OddRatioCalcAsm {
    private final String[] majorSNPReads;
    private final String[] minorSNPReads;
    private final double GeneExpressOddRatio;
    private final HashMap<Integer, String[]> aseSnpPos2BaseCount;

    /**
     * Constructor
     * @param majorSNPReads majorNc:majorCount:bamFileName:snpPos
     * @param minorSNPReads minorNc:minorCount:bamFileName:snpPos
     * @param aseSnpPos2BaseCount aseSnpPos -> [majorNc:majorCount:bam1, minorNc:minorCount:bam1,
     *                            majorNc:majorCount:bam2, minorNc:minorCount:bam2,...]
     */
    public OddRatioCalcAsm(String[] majorSNPReads, String[] minorSNPReads, double geneExpressOddRatio,
                           HashMap<Integer, String[]> aseSnpPos2BaseCount) {
        this.majorSNPReads = majorSNPReads;
        this.minorSNPReads = minorSNPReads;
        this.GeneExpressOddRatio = geneExpressOddRatio;
        this.aseSnpPos2BaseCount = aseSnpPos2BaseCount;
    }

    /**
     * calculate major allele LOR and the variance
     * @return {"LOR": [log orr ratios], "VAR": [variances]}
     */
    public HashMap<String, List<String>> getLogOddRatio() {
        List<String> logOR = new ArrayList<>();
        List<String> variance = new ArrayList<>();
        List<String> validSnpInfo = new ArrayList<>();
        String majorReads, minorReads;
        double logOddratio, var;
        for (int i = 0; i < this.majorSNPReads.length; i++) {
            majorReads = this.majorSNPReads[i];
            minorReads = this.minorSNPReads[i];
            logOddratio = this.calculateLogOddRatio(majorReads, minorReads);
            if (Double.isNaN(logOddratio)) {
                continue;
            }
            var = this.calculateVariance(majorReads, minorReads);
            logOR.add(String.valueOf(logOddratio));
            variance.add(String.valueOf(var));
            validSnpInfo.add(majorReads);
            validSnpInfo.add(minorReads);
        }
        HashMap<String, List<String>> calcRes = new HashMap<>();
        calcRes.put("LOR", logOR);
        calcRes.put("VAR", variance);
        calcRes.put("majorInfo", validSnpInfo);
    
        return calcRes;
    }

    /**
     * formula of LOR calculation
     *      yi = ln[y_major/(total-y_major)] - ln[0.5*total/(total-0.5*total)]
     *         = ln[y_major/(total-y_major)]
     *  majorAlleleReads List[majorNc:count:bam1:snp1, majorNc:count:bam1:snp1,
     *                       majorNc:count:bam2:snp1, majorNc:count:bam2:snp2,....]
     * @param majorAlleleReads majorNc:majorCount:bamFileName:snpPos
     * @param minorAlleleReads minorNc:minorCount:bamFileName:snpPos
     * @return LOR of a SNV site
     */
    private double calculateLogOddRatio(String majorAlleleReads, String minorAlleleReads) {
        double oddRatio = -1;
        try {
            String[] majorInfos = majorAlleleReads.split(":");
            String[] minorInfos = minorAlleleReads.split(":");
            int snpPos = Integer.parseInt(majorInfos[3]);
            String[] allele2Count = this.aseSnpPos2BaseCount.get(snpPos);
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
                oddRatio = (Double.parseDouble(majorInfos[1]) + 0.1) / (Double.parseDouble(minorInfos[1]) + 0.1) / this.GeneExpressOddRatio;
            } else if (aseMajor.equalsIgnoreCase(asmMinor) && aseMinor.equalsIgnoreCase(asmMajor)) {
                oddRatio = (Double.parseDouble(minorInfos[1]) + 0.1) / (Double.parseDouble(majorInfos[1]) + 0.1) / this.GeneExpressOddRatio;
            }
        } catch (Exception e) {
            return Double.NaN;
        }
        return Math.log(oddRatio);
    }

    /**
     * formula of LOR variance calculation
     *      s = 1/y_major + 1/(total-y_major) + 1/majorBackground + 1/minorBackground
     * @param majorAlleleReads MeRIP-seq INPUT data major allele reads count
     * @param minorAlleleReads MeRIP-seq INPUT data minor allele reads count
     * @return LOR variance of a SNV site
     */
    private double calculateVariance(String majorAlleleReads, String minorAlleleReads) {
        String[] majorInfos = majorAlleleReads.split(":");
        String[] minorInfos = minorAlleleReads.split(":");
        String bamFileName = majorInfos[2];
        int snpPos = Integer.parseInt(majorInfos[3]);
        String[] allele2Count = this.aseSnpPos2BaseCount.get(snpPos);
        String[] correspondAseSnpInfo = getCorrespondAseSnpInfo(allele2Count, bamFileName);
        assert correspondAseSnpInfo != null;
        return 1.0 / Integer.parseInt(majorInfos[1]) + 1.0 / Integer.parseInt(minorInfos[1]) + 1 / this.GeneExpressOddRatio + 1;
    }
    
    /**
     *
     * @param aseSnpInfo  [majorNc:majorCount:bam1, minorNc:minorCount:bam1, ...]
     * @param bamFileName bam1 || bam2...
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
}
