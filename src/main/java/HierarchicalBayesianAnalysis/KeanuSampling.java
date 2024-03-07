package HierarchicalBayesianAnalysis;

import io.improbable.keanu.Keanu;
import io.improbable.keanu.algorithms.NetworkSamples;
import io.improbable.keanu.algorithms.variational.optimizer.Optimizer;
import io.improbable.keanu.network.BayesianNetwork;
import io.improbable.keanu.network.KeanuProbabilisticModel;
import io.improbable.keanu.tensor.dbl.DoubleTensor;
import io.improbable.keanu.vertices.tensor.number.floating.dbl.probabilistic.GaussianVertex;
import io.improbable.keanu.vertices.tensor.number.floating.dbl.probabilistic.InverseGammaVertex;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author: tanglin
 * @date: 2022/07/18 17:35
 */
public class KeanuSampling {
    
    public static double sampleUseMaxPosterior(double df, double scaleParam, double[] variance, double[] oddRatio) {
        InverseGammaVertex tau = new InverseGammaVertex(df * 0.5, (double) df * scaleParam / 2);
        DoubleTensor varianceTensor = DoubleTensor.create(variance);
        double numerator = 0;
        double tauSample = tau.sample().scalar();
        for (int i = 0; i < varianceTensor.getLength(); i++) {
            numerator += oddRatio[i] / (variance[i] + Math.pow(tauSample, 2));
        }
        double muGene = numerator / varianceTensor.plus(tauSample).pow(-1).sum().scalar();
        double vGene = 1.0 / varianceTensor.plus(tauSample).pow(-1).sum().scalar();
        GaussianVertex u = new GaussianVertex(muGene, vGene);
        GaussianVertex[] theta = new GaussianVertex[oddRatio.length];
        GaussianVertex[] oddRatioVertex = new GaussianVertex[oddRatio.length];
        for (int i = 0; i < oddRatio.length; i++) {
            theta[i] = new GaussianVertex(u, tau);
            oddRatioVertex[i] = new GaussianVertex(theta[i], variance[i]);
        }
        for (int i = 0; i < oddRatioVertex.length; i++) {
            oddRatioVertex[i].observe(oddRatio[i]);
        }
        BayesianNetwork bayesNet = new BayesianNetwork(u.getConnectedGraph());
        Optimizer optimizer = Keanu.Optimizer.of(bayesNet);
        optimizer.maxAPosteriori();
        
        return u.getValue().scalar();
    }
    
    public static double[] sampleUseMh(double df, double scaleParam, double[] variance, double[] oddRatio, double initU, double initTau,
                                   double[] initTheta, int samplingTime, int burnIn) {
        InverseGammaVertex tau = new InverseGammaVertex(df * 0.5, df * scaleParam / 2);
        
        DoubleTensor varianceTensor = DoubleTensor.create(variance);
        DoubleTensor oddRatioTensor = DoubleTensor.create(oddRatio);
        
        DoubleTensor muGene = (oddRatioTensor.div(varianceTensor.plus(tau.pow(2).getValue())).sum()).div(varianceTensor.plus(tau.pow(2).getValue()).pow(-1).sum());
        DoubleTensor vGene = varianceTensor.plus(tau.pow(2).getValue()).pow(-1).sum();
        GaussianVertex u = new GaussianVertex(muGene, vGene.sqrt());
        
        DoubleTensor thetaV = (varianceTensor.pow(-1).plus(tau.pow(2).pow(-1).getValue())).pow(-1);
        DoubleTensor theta = (oddRatioTensor.div(varianceTensor).plus(muGene.div(tau.pow(2).getValue()))).times(thetaV);
        GaussianVertex oddRatioVertex = new GaussianVertex(theta, thetaV.sqrt());
        
        oddRatioVertex.observe(oddRatioTensor);
        
        u.setValue(initU);
        tau.setValue(initTau);
        oddRatioVertex.setValue(initTheta);
        
        KeanuProbabilisticModel model = new KeanuProbabilisticModel(u.getConnectedGraph());
        
        System.setProperty("io.improbable.keanu.util.status.StatusBar.disableStatusBar", "true");
        NetworkSamples posteriorSamples = Keanu.Sampling.MetropolisHastings.withDefaultConfig().getPosteriorSamples(
                model,
                model.getLatentVariables(),
                samplingTime
        ).drop(burnIn);
        return posteriorSamples.getDoubleTensorSamples(u).asTensor().asFlatDoubleArray();
    }
    
}
