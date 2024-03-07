package Threshold;

import Threshold.GradientDescent.GradientDescent;
import org.apache.commons.cli.*;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * @author: tanglin
 * @date: 2022/11/07 14:22
 * @Description:
 */
public class ThresholdComputeBySd {
    
    public static void main(String[] args) {
        Options options = new Options();
        CommandLine commandLine = null;
        try {
            commandLine = setCommandLine(args, options);
        } catch (ParseException pe) {
            System.err.println(pe.getMessage());
            System.exit(2);
        }
        String inputDirPath = new File(commandLine.getOptionValue("i")).getAbsolutePath();
        String outputFilePath = new File(commandLine.getOptionValue("o")).getAbsolutePath();
        int selectGeneNum = Integer.parseInt(commandLine.getOptionValue("n"));
        
        File dir = new File(inputDirPath);
        File[] distributionFiles = dir.listFiles();
        assert distributionFiles != null;
        String disName;
        List<RecordParamDTO> recordParamDTOS = new ArrayList<>();
        for (File distributionFile : distributionFiles) {
            disName = distributionFile.getName().split("\\.")[0];
            List<Double> scores = readScoreFromFile(distributionFile.getAbsolutePath(), selectGeneNum);
            RecordParamDTO recordParamDTO = getParam(scores, 0.05, disName);
            recordParamDTO.samplingMafs = scores;
            recordParamDTOS.add(recordParamDTO);
            System.out.println(disName + "successed!");
        }
        //saveParamRes(recordParamDTOS, new File(dir, "aseCategoryParam.txt").getAbsolutePath());
        saveParamRes(recordParamDTOS, outputFilePath);
    }
    
    private static void saveParamRes(List<RecordParamDTO> dtos, String resFilePath) {
        try(
                BufferedWriter bw = new BufferedWriter(new FileWriter(resFilePath))
        ){
            
            bw.write("distributionType" + "\t" + "gpdGamma(k)" + "\t" + "gpdSigma" + "\t" + "N" + "\t" + "Nt" + "\t" + "t" + "\t" + "samplingMafs");
            bw.newLine();
            for (RecordParamDTO dto : dtos) {
                StringBuilder sb = new StringBuilder();
                for (double samplingMaf : dto.samplingMafs) {
                    sb.append(String.format("%.5f", samplingMaf)).append(",");
                }
                String info = String.join("\t", dto.distributionType,
                        String.valueOf(dto.gpdGamma),
                        String.valueOf(dto.gpdSigma),
                        String.valueOf(dto.N),
                        String.valueOf(dto.Nt),
                        String.valueOf(dto.t),
                        sb.subSequence(0, sb.length() - 1));
                bw.write(info);
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private static RecordParamDTO getParam(List<Double> scores, double q, String distributionType) {
        int n = scores.size();
        Collections.sort(scores);
        double t = scores.get((int) ((1 - q - 0.01) * n) - 1);
        List<Double> gpdX = new ArrayList<>();
        for (Double score : scores) {
            double diff = score - t;
            if (diff > 0) {
                gpdX.add(diff);
            }
        }
        double[] X = new double[gpdX.size()];
        for (int i = 0; i < gpdX.size(); i++) {
            X[i] = gpdX.get(i);
        }
        GradientDescent gd = new GradientDescent(X);
        gd.optimize();
        
        RecordParamDTO dto = new RecordParamDTO();
        dto.distributionType = distributionType;
        dto.gpdGamma = gd.getGamma();
        dto.gpdSigma = gd.getSigma();
        dto.t = t;
        dto.N = n;
        dto.Nt = gpdX.size();
        return dto;
    }
 
    private static List<Double> readScoreFromFile(String filePath, int selectGeneNum) {
        List<Double> scores = new ArrayList<>();
        try (
                BufferedReader br = new BufferedReader(new FileReader(filePath))
        ){
            // 跳过表头
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                scores.add(Double.parseDouble(line.split("\t")[12]));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        Collections.shuffle(scores);
        return scores.subList(0, selectGeneNum);
    }
    
    private static CommandLine setCommandLine(String[] args, Options options) throws ParseException {
        Option option = new Option("o", "output", true, "result output file. Optional, default ./backgroundComputeRes.txt");
        option.setRequired(false);
        options.addOption(option);
        
        option = new Option("i", "input", true, "gene input file directory");
        option.setRequired(false);
        options.addOption(option);
    
        option = new Option("n", "selectNums", true, "the number of select genes");
        option.setRequired(false);
        options.addOption(option);
        
        CommandLineParser parser = new DefaultParser();
        
        return parser.parse(options, args);
    }
}
