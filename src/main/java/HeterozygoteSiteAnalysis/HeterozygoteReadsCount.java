package HeterozygoteSiteAnalysis;

import org.apache.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

public class HeterozygoteReadsCount {
    private final int readsCountThreshold;
    private final File peakCoveredSNPFile;
    
    private final HashMap<String, LinkedList<String>> peakMajorAlleleNc = new HashMap<>();
    private final Logger log;
    
    

    /**
     * Constructor
     * @param peakCoveredSnpFile file which record peak covered SNV information
     * @param readsCountThreshold reads count threshold
     * @param log log4j instance
     */
    public HeterozygoteReadsCount(String peakCoveredSnpFile, int readsCountThreshold, Logger log) {
        this.readsCountThreshold = readsCountThreshold;
        this.peakCoveredSNPFile = new File(peakCoveredSnpFile);
        this.log = log;
    }

    /**
     * allele reads count of each SNV sites under m6A signal range
     * isSnpBackground: Now whether to find the snp background number
     * @return [peak1: position1: [major: count: bamFile, minor: count:bamFile], position2:[major: count:bamFile, minor:count:bamFile]],
     * peak2:....]
     */
    public HashMap<String, HashMap<String, List<String>>> getMajorMinorHaplotype(boolean isSnpBackground) {
        HashMap<String, HashMap<String, List<String>>> majorMinorHaplotype = new HashMap<>();
        BufferedReader bfr = null;
        try {
            bfr = new BufferedReader(
                    new InputStreamReader(new FileInputStream(this.peakCoveredSNPFile))
            );
            String line = "";
            String[] info;
            String chr, position, peakStart, peakEnd, peakRange, geneId, majorNc, minorNc, bamFileName = null;
            int majorCount, minorCount;
            while (line != null) {
                line = bfr.readLine();
                if (line != null) {
                    if (line.startsWith("#"))
                        continue;
                    info = line.split("\t");
                    chr = info[0];
                    peakStart = info[1];
                    peakEnd = info[2];
                    position = info[3];
                    peakRange = peakStart +":"+peakEnd;
                    geneId = info[4];
                    majorNc = info[6];
                    minorNc = info[7];
                    majorCount = Integer.parseInt(info[8]);
                    minorCount = Integer.parseInt(info[9]);
                    if (!isSnpBackground) {
                        bamFileName = info[10];
                    }
                    
                    // the sum of two alleles >= 10
                    if (minorCount + majorCount < this.readsCountThreshold) {
                        continue;
                    }

                    //this.peakCoveredGene.put(chr+":"+peakStart+":"+peakEnd, geneId);
                    String label = String.join(":", chr, peakRange, geneId);
                    LinkedList<String> majorAlleleNcRecords = this.peakMajorAlleleNc.getOrDefault(label, new LinkedList<>());
                    String positionAndMajorNc = String.join(":", new String[]{position, majorNc});
                    /*
                        In the case of multiple BAM files, SNP loci under the same peak may appear in multiple BAM files,
                        resulting in duplicate recording of position and major nucleotide (majorNc).
                        Here, duplicate entries are removed to eliminate redundancy
                     */
                    if (!majorAlleleNcRecords.contains(positionAndMajorNc)) {
                        majorAlleleNcRecords.add(positionAndMajorNc);
                    }
                    this.peakMajorAlleleNc.put(label, majorAlleleNcRecords);
                    
                    HashMap<String, List<String>> peakSnp = majorMinorHaplotype.getOrDefault(label, new HashMap<>());
    
                    List<String> nucleotideReads = peakSnp.getOrDefault(position, new ArrayList<>());
                    // Code for processing wes data when multiple bam files are reserved
                    // it could need to change
                    String majorInfo = !isSnpBackground ? String.join(":", majorNc, String.valueOf(majorCount), bamFileName) :
                            String.join(":", majorNc, String.valueOf(majorCount));
                    String minorInfo = !isSnpBackground ? String.join(":", minorNc, String.valueOf(minorCount), bamFileName) :
                            String.join(":", minorNc, String.valueOf(minorCount));
                    nucleotideReads.add(majorInfo);
                    nucleotideReads.add(minorInfo);

                    peakSnp.put(position, nucleotideReads);
                    majorMinorHaplotype.put(label, peakSnp);
                }
            }
        } catch (IOException ie) {
            this.log.error("load file failed");
            this.log.error(ie.getMessage());
        } finally {
            if (bfr != null) {
                try {
                    bfr.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        return majorMinorHaplotype;
    }
}
