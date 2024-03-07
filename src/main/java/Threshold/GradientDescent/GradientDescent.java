package Threshold.GradientDescent;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author: tanglin
 * @date: 2022/11/07 10:24
 * @Description:
 */
public class GradientDescent {
    
    private double alpha = 0.03;
    private double diff = 0.0000001;
    private double[] x;
    private double bestTheta = 1;
    private int maxIterationNum = 5000;
    private double gamma;
    private double sigma;
    
    public GradientDescent() {
    }
    
    public GradientDescent(double[] x) {
        this.x = x;
        Arrays.sort(this.x);
    }
    
    public GradientDescent(double alpha, double diff, double[] x, int maxIterationNum, double bestTheta) {
        this.alpha = alpha;
        this.diff = diff;
        this.x = x;
        this.maxIterationNum = maxIterationNum;
        this.bestTheta = bestTheta;
        Arrays.sort(this.x);
    }
    
    public void optimize() {
        gradientDescent();
        double sum = 0;
        for (double s : x) {
            double v = 1 - bestTheta * s;
            if (v <= 0) {
                continue;
            }
            sum += Math.log(v);
        }
        gamma = -1.0 / x.length * sum;
        sigma = gamma / bestTheta;
    }
    
    public List<Integer> gdX = new ArrayList<>();
    public List<Double> gdY = new ArrayList<>();
    public void gradientDescent() {
        int iterationNum = 0;
        double fChange = bigGFunction(this.bestTheta);
        double fCurrent = bigGFunction(this.bestTheta);
        double temp;
        while (iterationNum < maxIterationNum && fChange > diff) {
            iterationNum ++;
            if (iterationNum % 10 == 0) {
                System.out.println("have computed" + iterationNum + "次");
            }
            bestTheta = this.bestTheta - alpha * bigGFunctionDerivative(bestTheta);
            temp = bigGFunction(bestTheta);
            fChange = Math.abs(fCurrent - temp);
            fCurrent = temp;
            gdX.add(iterationNum);
            gdY.add(temp);
        }
        System.out.println("diff is " + fChange);
    }
    
    private double bigGFunction(double theta) {
        int n = x.length;
        
        double gThetaSum = calculateGThetaSum(theta);
        
        double sum = 0;
        for (int i = 0; i < n; i++) {
            double gFunctionValue = gFunction(gThetaSum, theta, i);
            if (gFunctionValue == 0) {
                sum += 0;
            } else {
                sum += (2 * i - 1) * Math.log(gFunction(gThetaSum, theta, i)) + (2 * n + 1 - 2 * i) * Math.log(1 - gFunction(gThetaSum, theta, i));
            }
        }
        return -n - 1.0 / n * sum;
    }
    
    private double bigGFunctionDerivative(double theta) {
        double gThetaSum = calculateGThetaSum(theta);
        
        double sum = 0;
        int n = x.length;
        for (int i = 0; i < n; i++) {
            double gFunctionValue = gFunction(gThetaSum, theta, i);
            if (gFunctionValue == 0) {
                sum += 0;
            } else {
                sum += (2 * i - 1) * gFunctionDerivative(theta, i) / gFunction(gThetaSum, theta, i) -
                        (2 * n + 1 - 2 * i) * gFunctionDerivative(theta, i) / (1 - gFunction(gThetaSum, theta, i));
            }
        }
        return -1.0 / n * sum;
    }
    
    private double gFunction(double gThetaSum, double theta, int index) {
        double exponent = -x.length / gThetaSum;
        double temp = 1- theta * x[index];
        if (temp <= 0) {
            return 0;
        }
        return 1 - Math.pow((1 - theta * x[index]), exponent);
    }
    
    private double gFunctionDerivative(double theta, int index) {
        if (1 - theta * x[index] <= 0) {
            return  0;
        }
        double sum1 = 0;
        double sum2 = 0;
        int n = x.length;
        for (double s : x) {
            double temp = 1 - theta * s;
            if (temp <= 0) {
                continue;
            }
            sum1 += Math.log(1 - theta * s);
            sum2 -= s / (1 - theta * s);
        }
        return -1 * Math.pow((1 - theta * x[index]), -n / sum1) *
                (n * x[index] / ((1 - theta * x[index]) * sum1) + n * Math.log(1 - theta * x[index]) * sum2 / Math.pow(sum1, 2));
    }
    
    private double calculateGThetaSum(double theta) {
        double gThetaSum = 0;
        for (double t : x) {
            if (1 - theta * t <= 0) {
                gThetaSum += 0;
            } else {
                gThetaSum += Math.log(1 - theta * t);
            }
        }
        return gThetaSum;
    }
    
    public double getGamma() {
        return gamma;
    }
    
    public double getSigma() {
        return sigma;
    }
}
