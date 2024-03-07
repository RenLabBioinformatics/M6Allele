package TreeUtil;

import htsjdk.tribble.index.interval.Interval;

/**
 * @author: tanglin
 * @date: 2022/10/26 16:58
 * @Description: gtf interval tree、exon interval tree and peak interval tree node
 */
public class IntervalTreeNode extends Interval {
    public String geneName, geneId, strand;
    public int start, end;
    
    public IntervalTreeNode(int start, int end) {
        super(start, end);
        this.start = start;
        this.end = end;
    }
    
    public IntervalTreeNode(int start, int end, String geneId) {
        super(start, end);
        this.geneId = geneId;
        this.start = start;
        this.end = end;
    }
    
    public IntervalTreeNode(int start, int end, String geneId, String geneName, String strand) {
        super(start, end);
        this.start = start;
        this.end = end;
        this.geneId = geneId;
        this.geneName = geneName;
        this.strand = strand;
    }
}
