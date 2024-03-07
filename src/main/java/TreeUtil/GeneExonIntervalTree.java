package TreeUtil;

import htsjdk.tribble.index.interval.IntervalTree;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

public class GeneExonIntervalTree {
    private final String gtfFile;
    private HashMap<String, HashMap<String, IntervalTree>> geneExonIntervalTree;

    public GeneExonIntervalTree(String gtfFile) {
        this.gtfFile = gtfFile;
    }

    public void generateExonTree() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.gtfFile)));
            this.geneExonIntervalTree = new HashMap<>();
            String line = "", chrNum, geneId, geneName, strand;
            String[] info, geneInfo;
            int exonStart, exonEnd;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    if (!info[2].equalsIgnoreCase("exon"))
                        continue;
                    chrNum = info[0];
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];
                    exonStart = Integer.parseInt(info[3]);
                    exonEnd = Integer.parseInt(info[4]);
                    strand = info[6];
                    IntervalTreeNode gitn = new IntervalTreeNode(exonStart, exonEnd, geneId, geneName, strand);

                    HashMap<String, IntervalTree> gene2ExonTree = this.geneExonIntervalTree.getOrDefault(chrNum, new HashMap<>());
                    IntervalTree it = gene2ExonTree.getOrDefault(geneId, new IntervalTree());
                    it.insert(gitn);
                    gene2ExonTree.put(geneId, it);
                    this.geneExonIntervalTree.put(chrNum, gene2ExonTree);
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

    /**
     * get gene name
     * @param recordInfo GTF information
     * @return gene name
     */
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

    public HashMap<String, HashMap<String, IntervalTree>> getGeneExonIntervalTree() {
        return this.geneExonIntervalTree;
    }
}
