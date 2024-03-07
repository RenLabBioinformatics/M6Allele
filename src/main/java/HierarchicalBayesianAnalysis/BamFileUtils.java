package HierarchicalBayesianAnalysis;

import htsjdk.samtools.*;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author: tanglin
 * @date: 2022/04/21 20:21
 * @Description:
 */
public class BamFileUtils {
    /**
     * Find the base and the corresponding count for the major and minor of the specified snp in the bam file
     */
    public static String getBase2CountInBam(String bam, String bai, String chr, int snpPosition) {
        SamInputResource samInputResource = SamInputResource.of(bam).index(new File(bai));
        try (SamReader sr = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS).
                validationStringency(ValidationStringency.SILENT).open(samInputResource)) {
            Map<String, Integer> base2CountMap = new HashMap<>();
            int sequenceIndex = sr.getFileHeader().getSequenceDictionary().getSequenceIndex(chr);
            QueryInterval[] queryParam = new QueryInterval[] {new QueryInterval(sequenceIndex, snpPosition, snpPosition)};
            
            SAMRecordIterator samIterator = sr.queryOverlapping(queryParam);
            SAMRecord samRecord;
            while (samIterator.hasNext()) {
                samRecord = samIterator.next();
                List<AlignmentBlock> alignmentBlocks = samRecord.getAlignmentBlocks();
                String read = samRecord.getReadString();
                int index = -1;
                int span = 0;
                for (int i = 0; i < alignmentBlocks.size(); i++) {
                    int blockReferenceStart = alignmentBlocks.get(i).getReferenceStart();
                    int blockLength = alignmentBlocks.get(i).getLength();
                    int blockReferenceEnd = blockReferenceStart + blockLength - 1;
                    if (blockReferenceStart <= snpPosition && snpPosition <= blockReferenceEnd) {
                        span = snpPosition - blockReferenceStart;
                        index = i;
                        break;
                    }
                }
                if (index != -1) {
                    int readStart = alignmentBlocks.get(index).getReadStart();
                    String base = read.charAt(readStart + span - 1) + "";
                    Integer count = base2CountMap.getOrDefault(base, 0);
                    count++;
                    base2CountMap.put(base, count);
                }
            }
            samIterator.close();
            List<Map.Entry<String, Integer>> baseCount = base2CountMap.entrySet().stream().sorted(Map.Entry.<String, Integer>comparingByValue().reversed())
                    .collect(Collectors.toList());
            if (baseCount.size() < 2) {
                return null;
            }
            String majorNc = baseCount.get(0).getKey();
            int majorAlleleCount = baseCount.get(0).getValue();
            String minorNc = baseCount.get(1).getKey();
            int minorAlleleCount = baseCount.get(1).getValue();
            return String.join(":", majorNc, String.valueOf(majorAlleleCount), minorNc, String.valueOf(minorAlleleCount));
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }
    
    public static void close(SamReader sr) {
        if (sr != null) {
            try {
                sr.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
