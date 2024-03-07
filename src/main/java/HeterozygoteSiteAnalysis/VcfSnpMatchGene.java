package HeterozygoteSiteAnalysis;

import TreeUtil.GTFIntervalTree;
import TreeUtil.GeneExonIntervalTree;
import TreeUtil.IntervalTreeNode;
import HierarchicalBayesianAnalysis.BamFileUtils;
import htsjdk.tribble.index.interval.Interval;
import htsjdk.tribble.index.interval.IntervalTree;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 * match VCF file mutation record to corresponding gene
 */
public class VcfSnpMatchGene {
    protected String vcfFile;
    private final int readsCoverageThreshold;
    // GTF interval tree for each chromosome
    private HashMap<String, htsjdk.tribble.index.interval.IntervalTree> gtfIntervalTree;
    private HashMap<String, HashMap<String, htsjdk.tribble.index.interval.IntervalTree>> geneExonIntervalTree;
    // geneAlleleReads = {"geneId->geneName": {pos1: [majorAllele:count:bam1, minorAllele: count1:bam1]}, ...}
    private final HashMap<String, HashMap<Integer, String[]>> geneAlleleReads = new HashMap<>();
    public VcfSnpMatchGene(String vcfFile, String gtfFile, int readsCoverageThreshold) {
        this.vcfFile = vcfFile;
        Logger logger = Logger.getLogger(VcfSnpMatchGene.class);
        logger.info("start processing gtf file");
        GTFIntervalTree git = new GTFIntervalTree(gtfFile);
        git.parseGTFFile();
        GeneExonIntervalTree geit = new GeneExonIntervalTree(gtfFile);
        geit.generateExonTree();
        this.gtfIntervalTree = git.getGtfIntervalTrees();
        this.geneExonIntervalTree = geit.getGeneExonIntervalTree();
        this.readsCoverageThreshold = readsCoverageThreshold;
    }

    public HashMap<String, HashMap<Integer, String[]>> getGeneAlleleReads() {
        return this.geneAlleleReads;
    }
    
    
    public void parseBamFile(List<String> bamFile, List<String> baiFile) {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.vcfFile)));
            String line = "", chrNum;
            int position;
            IntervalTree it;
            List<Interval> itn;
            HashMap<String, IntervalTree> eit;
            IntervalTreeNode gitn;
            String[] info;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    
                    List<String> majorNcs = new ArrayList<>();
                    List<String> minorNcs = new ArrayList<>();
                    
                    List<Integer> majorAlleleCounts = new ArrayList<>();
                    List<Integer> minorAlleleCounts = new ArrayList<>();
                    
                    List<String> bamRecords = new ArrayList<>();
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    position = Integer.parseInt(info[1]);
                    it = this.gtfIntervalTree.getOrDefault(chrNum, null);
                    if (it == null)
                        continue;
                    eit = this.geneExonIntervalTree.getOrDefault(chrNum, null);
                    for (int i = 0; i < bamFile.size(); i++) {
                        String base2CountInBam = BamFileUtils.getBase2CountInBam(bamFile.get(i), baiFile.get(i), chrNum, position);
                        if (base2CountInBam == null) {
                            continue;
                        }
                        String[] baseCount = base2CountInBam.split(":");
                        // each allele >= 2  the sum of two allele >= 10
                        if (Integer.parseInt(baseCount[3]) + Integer.parseInt(baseCount[1]) < this.readsCoverageThreshold ||
                            Integer.parseInt(baseCount[3]) < 2 || Integer.parseInt(baseCount[1]) < 2) {
                            continue;
                        }
                        majorNcs.add(baseCount[0]);
                        majorAlleleCounts.add(Integer.parseInt(baseCount[1]));
                        minorNcs.add(baseCount[2]);
                        minorAlleleCounts.add(Integer.parseInt(baseCount[3]));
                        bamRecords.add("bam" + (i + 1));
                    }
                    // locate the SNV site on particular gene using GTF interval tree
                    itn = it.findOverlapping(new Interval(position, position));
                    if (CollectionUtils.isEmpty(itn) || CollectionUtils.isEmpty(majorNcs) || CollectionUtils.isEmpty(majorAlleleCounts) ||
                        CollectionUtils.isEmpty(minorNcs) || CollectionUtils.isEmpty(minorAlleleCounts)) {
                        continue;
                    }
                    for (Interval interval : itn) {
                        gitn = (IntervalTreeNode) interval;
                        this.recursiveSearch(gitn, eit, majorNcs, minorNcs, majorAlleleCounts, minorAlleleCounts, position, bamRecords);
                    }
                }
            }
        } catch (IOException ie) {
            ie.printStackTrace();
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            this.gtfIntervalTree = null;
            this.geneExonIntervalTree = null;
        }
    }

    /**
     * renew geneAlleleReads
     * @param geneId gene id
     * @param geneName gene name
     * @param majorAllele major allele nucleotide
     * @param minorAllele minor allele nucleotide
     * @param mutPosition mutation position
     * @param majorAlleleCount major allele reads count
     * @param minorAlleleCount minor allele reads count
     */
    private void renewGeneAlleleReads(String geneId, String geneName, List<String> majorAllele, List<String> minorAllele,
                                      int mutPosition, List<Integer> majorAlleleCount,
                                      List<Integer> minorAlleleCount, List<String> bamRecords) {
        String label = String.join("->", new String[]{geneId, geneName});
        HashMap<Integer, String[]> geneSnp = this.geneAlleleReads.getOrDefault(label, new HashMap<>());
        List<String> referencesAndAlternatives = new ArrayList<>();
        for (int i = 0; i < majorAllele.size(); i++) {
            String reference = String.join(":", new String[]{majorAllele.get(i), String.valueOf(majorAlleleCount.get(i)), bamRecords.get(i)});
            String alternative = String.join(":", new String[]{minorAllele.get(i), String.valueOf(minorAlleleCount.get(i)), bamRecords.get(i)});
            referencesAndAlternatives.add(reference);
            referencesAndAlternatives.add(alternative);
        }
        geneSnp.put(mutPosition, referencesAndAlternatives.toArray(new String[0]));
        this.geneAlleleReads.put(label, geneSnp);
    }

    /**
     * recursive search
     * @param gitn GTFIntervalTreeNode instance on GTFIntervalTree
     * @param exonIntervalTree IntervalTree records genes' exon region
     * @param majorAllele major allele nucleotide
     * @param minorAllele minor allele nucleotide
     * @param majorAlleleCount major allele reads count
     * @param minorAlleleCount minor allele reads count
     * @param position mutation position
     */
    public void recursiveSearch(IntervalTreeNode gitn, HashMap<String, IntervalTree> exonIntervalTree,
                                List<String> majorAllele, List<String> minorAllele,
                                List<Integer> majorAlleleCount, List<Integer> minorAlleleCount, int position, List<String> bamRecords) {
        String geneName = gitn.geneName;
        String geneId = gitn.geneId;
        
        IntervalTree intervalTreeExon = exonIntervalTree.getOrDefault(geneId, null);
        if (intervalTreeExon != null) {
            
            List<Interval> potentialExon = intervalTreeExon.findOverlapping(new Interval(position, position));
            if (CollectionUtils.isNotEmpty(potentialExon)) {
                this.renewGeneAlleleReads(geneId, geneName, majorAllele, minorAllele, position, majorAlleleCount, minorAlleleCount, bamRecords);
            }
        }
    }
}
