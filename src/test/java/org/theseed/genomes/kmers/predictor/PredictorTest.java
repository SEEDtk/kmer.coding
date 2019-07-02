/**
 *
 */
package org.theseed.genomes.kmers.predictor;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import org.theseed.genomes.Contig;
import org.theseed.genomes.Genome;
import org.theseed.genomes.kmers.DnaKmer;
import org.theseed.genomes.kmers.SequenceDnaKmers;
import org.theseed.locations.Frame;
import org.theseed.locations.LocationList;
import org.theseed.utils.CountMap;

import com.github.cliftonlabs.json_simple.JsonException;

import junit.framework.TestCase;

/**
 * @author Bruce Parrello
 *
 */

public class PredictorTest extends TestCase {

    /**
     * @param name
     */
    public PredictorTest(String name) {
        super(name);
    }

    /**
     * Tests for FramePredictor
     * @throws IOException
     */
    public void testFramePredictor() throws IOException {
        FramePredictor testPred = new FramePredictor("kmers.tbl");
        DnaKmer testKmer = new DnaKmer("gacgggcgtgtagac");
        assertEquals("Test kmer in wrong frame.", Frame.M0, testPred.frameOf(testKmer));
        assertEquals("Non-coding kmer in wrong frame.", Frame.F0, testPred.frameOf("gacgggcggtgtgtg"));
        assertEquals("Plus-one kmer in wrong frame.", Frame.P1, testPred.frameOf("gacgggctacacatt"));
    }

    public void testFramePredictions() throws NumberFormatException, IOException, JsonException {
        // We will track counts in here.
        CountMap<DnaKmer> testCounts = new CountMap<DnaKmer>();
        // Load the predictor.
        FramePredictor testPred = new FramePredictor("\\Users\\Bruce\\Documents\\FIG\\SEEDtk\\Data\\TestKmers\\kmers.tbl");
        String[] genomes = new String[] { "bin3.gto", "bin4.gto" };
        for (String genome : genomes) {
            // Load the genome.
            Genome myGto = new Genome(genome);
            Map<String, LocationList> gtoMap = LocationList.createGenomeCodingMap(myGto);
            // Loop through the contigs.
            Collection<Contig> allContigs = myGto.getContigs();
            for (Contig contig : allContigs) {
                SequenceDnaKmers contigKmers = new SequenceDnaKmers(contig.getSequence());
                LocationList contigLocs = gtoMap.get(contig.getId());
                while (contigKmers.nextKmer()) {
                    int pos = contigKmers.getPos();
                    Frame predicted = testPred.frameOf(contigKmers);
                    if (predicted != Frame.XX) {
                        Frame kmerFrame = contigLocs.computeRegionFrame(pos, pos + DnaKmer.getSize() - 1);
                        if (kmerFrame != Frame.XX) {
                            // Get an immutable copy of the kmer to use as a key.
                            DnaKmer kmer = contigKmers.getCopy();
                            if (! kmerFrame.equals(predicted)) {
                                testCounts.setBad(kmer);
                            } else {
                                testCounts.setGood(kmer);
                            }
                        }
                    }
                }
            }
        }
        // Loop through the results, looking for problems.
        Collection<DnaKmer> results = testCounts.keys();
        CountMap<Frame> frameCounts = new CountMap<Frame>();
        int badHits = 0;
        int goodHits = 0;
        for (DnaKmer kmer : results) {
            int good = testCounts.good(kmer);
            int bad = testCounts.bad(kmer);
            badHits += bad;
            goodHits += good;
            if (bad > 0) {
                frameCounts.setBad(testPred.frameOf(kmer));
            } else {
                frameCounts.setGood(testPred.frameOf(kmer));
            }
        }
        System.err.format("Total hits = %d good, %d bad.%n", goodHits, badHits);
        System.err.println("Summary of kmers by frame.");
        for (Frame frm : Frame.all) {
            System.err.format("%-4s %8d %8d%n", frm, frameCounts.good(frm), frameCounts.bad(frm));
        }
    }

}
