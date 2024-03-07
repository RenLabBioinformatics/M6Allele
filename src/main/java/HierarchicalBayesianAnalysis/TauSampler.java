package HierarchicalBayesianAnalysis;

public class TauSampler extends MHSampling {
    private final ScaledInvChiSquareDistribution priorTau;
    
    /**
     * Constructor
     * @param df number of chi-squared degrees of freedom
     */
    public TauSampler(double df, double scaleParam) {
        super();
        this.priorTau = new ScaledInvChiSquareDistribution(df, scaleParam);
    }

    /**
     * sampling a new tau from its priority distribution
     * @return randomly sample tau
     */
    public double randomInit() {
        return this.priorTau.sample();
    }
    
}
