package HeterozygoteSiteAnalysis;

import TreeUtil.GTFIntervalTree;
import TreeUtil.GeneExonIntervalTree;
import TreeUtil.IntervalTreeNode;
import HierarchicalBayesianAnalysis.BamFileUtils;
import htsjdk.tribble.index.interval.Interval;
import htsjdk.tribble.index.interval.IntervalTree;
import org.apache.commons.collections4.CollectionUtils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class GeneSNVRecord {
    private final String gtfFile;
    private final String vcfFile;
    private final String outputFile;
    private HashMap<String, IntervalTree> geneIntervalTree;
    private HashMap<String, HashMap<String, IntervalTree>> geneExonIntervalTree;
    

    public GeneSNVRecord(String gtfFile, String vcfFile, String outputFile) {
        this.gtfFile = gtfFile;
        this.vcfFile = vcfFile;
        this.outputFile = outputFile;
    }
    
    
    public void locateSnv(List<String> bamFile, List<String> baiFile) {
        this.parseGTFFile();
        this.parseVCFFile(bamFile, baiFile);
    }

    private void parseGTFFile() {
        GTFIntervalTree git = new GTFIntervalTree(this.gtfFile);
        git.parseGTFFile();
        this.geneIntervalTree = git.getGtfIntervalTrees();
        GeneExonIntervalTree geit = new GeneExonIntervalTree(this.gtfFile);
        geit.generateExonTree();
        this.geneExonIntervalTree = geit.getGeneExonIntervalTree();
    }
    
    private void parseVCFFile(List<String> bamFile, List<String> baiFile) {
        BufferedReader bfr = null;
        BufferedWriter bfw = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.vcfFile)));
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(this.outputFile)));
            bfw.write(String.join("\t", new String[] {"chr", "geneId", "geneName", "position",
                    "majorAllele", "minorAllele",
                    "majorAlleleCount", "minorAlleleCount"}));
            bfw.newLine();
            String readIn = "", writeOut, chrNum, geneId, geneName, majorNc, minorNc;
            int position, majorAlleleCount, minorAlleleCount;
            String[] info;
            List<IntervalTreeNode> potentialLocateGenes;
            while (readIn != null) {
                readIn = bfr.readLine();
                if (readIn != null) {
                    if (readIn.startsWith("#"))
                        continue;
                    info = readIn.split("\t");
                    chrNum = info[0];
                    IntervalTree it = this.geneIntervalTree.getOrDefault(chrNum, null);
                    if (it == null)
                        continue;
                    position = Integer.parseInt(info[1]);
                    
                    List<Interval> selectGenes = it.findOverlapping(new Interval(position, position));
                    if (CollectionUtils.isEmpty(selectGenes)) {
                        continue;
                    }
                    
                    potentialLocateGenes = this.searchLocateGene(selectGenes, chrNum, position);
                    for (int i = 0; i < baiFile.size(); i++) {
                        String baseToCount = BamFileUtils.getBase2CountInBam(bamFile.get(i), baiFile.get(i), chrNum, position);
                        if (baseToCount == null) {
                            continue;
                        }
                        String[] baseCount = baseToCount.split(":");
                        majorNc = baseCount[0];
                        majorAlleleCount = Integer.parseInt(baseCount[1]);
                        minorNc = baseCount[2];
                        minorAlleleCount = Integer.parseInt(baseCount[3]);
                        
                        for (IntervalTreeNode node: potentialLocateGenes) {
                            geneId = node.geneId;
                            geneName = node.geneName;
                            writeOut = String.join("\t", new String[]{chrNum, geneId, geneName, info[1], majorNc, minorNc,
                                    String.valueOf(majorAlleleCount), String.valueOf(minorAlleleCount)});
                            bfw.write(writeOut);
                            bfw.newLine();
                        }
                    }
                }
            }
        } catch (IOException ie) {
            ie.printStackTrace();
            System.exit(2);
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
            if (bfw != null) {
                try {
                    bfw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private List<IntervalTreeNode> searchLocateGene(List<Interval> selectGenes, String chrNum, int position) {
        List<IntervalTreeNode> potentialGenes = new ArrayList<>();
        for (Interval interval : selectGenes) {
            IntervalTreeNode gene = (IntervalTreeNode) interval;
            String geneId = gene.geneId;
            boolean inExon = this.ifInExon(chrNum, geneId, position);
            if (inExon) {
                potentialGenes.add(gene);
            }
        }
        return potentialGenes;
    }

    private boolean ifInExon(String chrNum, String geneId, int position) {
        HashMap<String, IntervalTree> exonTreeMap = this.geneExonIntervalTree.getOrDefault(chrNum, null);
        if (exonTreeMap == null) {
            return false;
        }
        IntervalTree it = exonTreeMap.getOrDefault(geneId, null);
        if (it == null) {
            return false;
        }
        List<Interval> res = it.findOverlapping(new Interval(position, position));
        return CollectionUtils.isNotEmpty(res);
    }

}
