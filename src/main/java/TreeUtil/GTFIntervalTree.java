package TreeUtil;

import htsjdk.tribble.index.interval.IntervalTree;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

/**
 * @author: tanglin
 * @date: 2022/10/26 17:02
 * @Description:
 */
public class GTFIntervalTree {
    private final String gtfFile;
    private final HashMap<String, IntervalTree> gtfIntervalTrees = new HashMap<>();
    
    public GTFIntervalTree(String gtfFile) {
        this.gtfFile = gtfFile;
    }
    
    /**
     * parse GTF file, build interval tree for each chromosome
     */
    public void parseGTFFile() {
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(new InputStreamReader(new FileInputStream(this.gtfFile)));
            String line = "", chrNum, geneName, geneId, componentType, strand;
            String[] info, geneInfo;
            int start, end;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chrNum = info[0];
                    componentType = info[2];
                    if (!componentType.equals("gene"))
                        continue;
                    start = Integer.parseInt(info[3]);
                    end = Integer.parseInt(info[4]);
                    strand = info[6];
                    geneInfo = this.getGeneInfo(info[8]);
                    geneId = geneInfo[0];
                    geneName = geneInfo[1];
                    
                    IntervalTreeNode gtfNode = new IntervalTreeNode(start, end, geneId, geneName, strand);
                    IntervalTree it = this.gtfIntervalTrees.getOrDefault(chrNum, new IntervalTree());
                    it.insert(gtfNode);
                    this.gtfIntervalTrees.put(chrNum, it);
                }
            }
            bfr.close();
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
    
    /**
     * return interval tree
     * @return interval tree
     */
    public HashMap<String, IntervalTree> getGtfIntervalTrees() {
        return this.gtfIntervalTrees;
    }
}
