package HeterozygoteSiteAnalysis;

import TreeUtil.GTFIntervalTree;
import TreeUtil.GeneExonIntervalTree;
import TreeUtil.IntervalTreeNode;
import HierarchicalBayesianAnalysis.BamFileUtils;
import htsjdk.tribble.index.interval.Interval;
import htsjdk.tribble.index.interval.IntervalTree;
import org.apache.commons.collections4.CollectionUtils;

import java.io.*;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.List;

public class PeakWithSNV {
    private final String gtfFile;
    private final String bedFile;
    private final String vcfFile;
    private final String outputFile;
    private HashMap<String, htsjdk.tribble.index.interval.IntervalTree> geneIntervalTree;
    // chr -> geneId -> peakTree;  chr -> geneId -> exonIntervalTree
    private HashMap<String, HashMap<String, htsjdk.tribble.index.interval.IntervalTree>> peakTreeMap, geneExonIntervalTree;
    private HashMap<String, String> geneNames;

    public PeakWithSNV(String gtfFile, String bedFile, String vcfFile, String outputFile) {
        this.gtfFile = gtfFile;
        this.bedFile = bedFile;
        this.vcfFile = vcfFile;
        this.outputFile = outputFile;
    }
    
    
    public void locateSnvInPeak(List<String> bamFile, List<String> baiFile) {
        this.buildGTFIntervalTree();
        this.buildPeakIntervalTreeWithoutSplit();
        this.parseGTFFile();
        this.parseVCFFile(bamFile, baiFile);
    }

    private void buildGTFIntervalTree() {
        GTFIntervalTree git = new GTFIntervalTree(this.gtfFile);
        git.parseGTFFile();
        this.geneIntervalTree = git.getGtfIntervalTrees();
    }
    
    private void buildPeakIntervalTreeWithoutSplit() {
        BufferedReader bfr = null;
        String chrNum, geneId;
        int start, end, peakCenter;
        List<Interval> potentialLocateGenes;
        this.peakTreeMap = new HashMap<>();
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.bedFile)));
            String line = "";
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#")) {
                        continue;
                    }
                    info = line.split("\t");
                    // BED format file must contains chr, peakStart, peakEnd
                    chrNum = info[0];
                    start = Integer.parseInt(new BigDecimal(info[1]).toPlainString());
                    end = Integer.parseInt(new BigDecimal(info[2]).toPlainString());
                    peakCenter = (int) Math.round((start + end) * 0.5);
                    
                    IntervalTree it = this.geneIntervalTree.getOrDefault(chrNum, null);
                    if (it == null)
                        break;
                    potentialLocateGenes = it.findOverlapping(new Interval(peakCenter, peakCenter));
                    if (CollectionUtils.isEmpty(potentialLocateGenes)) {
                        continue;
                    }
                    
                    // chr -> geneId -> peakTree
                    HashMap<String, IntervalTree> chrTree = peakTreeMap.getOrDefault(chrNum, new HashMap<>());
                    
                    for (Interval node: potentialLocateGenes) {
                        IntervalTreeNode potentialGene = (IntervalTreeNode) node;
                        geneId = potentialGene.geneId;
                        IntervalTreeNode peakNode = new IntervalTreeNode(start, end);
                        IntervalTree tree = chrTree.getOrDefault(geneId, new IntervalTree());
                        tree.insert(peakNode);
                        chrTree.put(geneId, tree);
                    }
                    peakTreeMap.put(chrNum, chrTree);
                    
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
        }
    }
    

    private void parseGTFFile() {
        this.geneIdMatchName();
        this.generateExonTree();
    }

    private void geneIdMatchName() {
        this.geneNames = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.gtfFile)));
            String line = "", geneId, geneName;
            String[] info, geneInfo;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    if (!info[2].equals("gene"))
                        continue;
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];
                    this.geneNames.put(geneId, geneName);
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
        }
    }

    private void generateExonTree() {
        GeneExonIntervalTree geit = new GeneExonIntervalTree(this.gtfFile);
        geit.generateExonTree();
        this.geneExonIntervalTree = geit.getGeneExonIntervalTree();
    }

    private String[] getGeneInfo(String recordInfo) {
        String[] info = recordInfo.split("; ");
        String geneName = "unknown", geneId = null;
        for (String s: info) {
            if (s.startsWith("gene_id")) {
                String[] name = s.split(" ");
                geneId = name[1].substring(1, name[1].length() -1);
            }
            if (s.startsWith("gene_name")) {
                String[] name = s.split(" ");
                geneName = name[1].substring(1, name[1].length() -1);
            }
        }

        return new String[] {geneId, geneName};
    }

    
    private void parseVCFFile(List<String> bamFile, List<String> baiFile) {
        BufferedReader bfr = null;
        BufferedWriter bfw = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.vcfFile)));
            bfw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(this.outputFile)));
            bfw.write(String.join("\t", new String[] {"#chr", "peakStart", "peakEnd", "mutatePosition",
                    "geneId", "geneName", "majorAllele", "minorAllele",
                    "majorAlleleCount", "minorAlleleCount", "bamFile"}));
            bfw.newLine();
            String readIn = "", writeOut, chrNum, majorNc, minorNc;
            int position, majorAlleleCount, minorAlleleCount;
            String[] info;
            IntervalTree geneTree, peakTree;
            while (readIn != null) {
                readIn = bfr.readLine();
                if (readIn != null) {
                    if (readIn.startsWith("#"))
                        continue;
                    info = readIn.split("\t");
                    chrNum = info[0];
                    position = Integer.parseInt(info[1]);
                    // gene -> peakTree
                    HashMap<String, IntervalTree> trees = this.peakTreeMap.getOrDefault(chrNum, null);
                    if (trees == null) {
                        continue;
                    }
                    geneTree = this.geneIntervalTree.getOrDefault(chrNum, null);
                    if (geneTree == null) {
                        continue;
                    }
                    for (int i = 0; i < bamFile.size(); i++) {
                        String bamFileName = "bam" + (i + 1);
                        String base2CountInBam = BamFileUtils.getBase2CountInBam(bamFile.get(i), baiFile.get(i), chrNum, position);
                        if (base2CountInBam == null) {
                            continue;
                        }
                        String[] baseCount = base2CountInBam.split(":");
                        majorNc = baseCount[0];
                        majorAlleleCount = Integer.parseInt(baseCount[1]);
                        minorNc = baseCount[2];
                        minorAlleleCount = Integer.parseInt(baseCount[3]);
                        
                        List<Interval> potentialGenes = geneTree.findOverlapping(new Interval(position, position));
                        for (Interval node : potentialGenes) {
                            IntervalTreeNode gene = (IntervalTreeNode) node;
                            
                            if (!this.ifInExon(chrNum, gene.geneId, position)) {
                                continue;
                            }
                            
                            peakTree = trees.getOrDefault(gene.geneId, null);
                            if (peakTree == null) {
                                continue;
                            }
                            
                            List<Interval> potentialPeaks = peakTree.findOverlapping(new Interval(position, position));
                            if (CollectionUtils.isNotEmpty(potentialPeaks)) {
                                for (Interval peakNode : potentialPeaks) {
                                    IntervalTreeNode peak = (IntervalTreeNode) peakNode;
                                    writeOut = String.join("\t",
                                            new String[] {chrNum, String.valueOf(peak.start), String.valueOf(peak.end),
                                                    info[1], gene.geneId, gene.geneName, majorNc, minorNc, String.valueOf(majorAlleleCount),
                                                    String.valueOf(minorAlleleCount), bamFileName});
                                    bfw.write(writeOut);
                                    bfw.newLine();
                                }
                            }
                        }
                    }
                }
            }
            this.geneNames = null;
            this.peakTreeMap = null;
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
    
    private boolean ifInExon(String chrNum, String geneId, int position) {
        HashMap<String, IntervalTree> treeMap = this.geneExonIntervalTree.getOrDefault(chrNum, null);
        if (treeMap == null)
            return false;
        IntervalTree it = treeMap.getOrDefault(geneId, null);
        if (it == null)
            return false;
        List<Interval> itn = it.findOverlapping(new Interval(position, position));

        return CollectionUtils.isNotEmpty(itn);
    }
}
