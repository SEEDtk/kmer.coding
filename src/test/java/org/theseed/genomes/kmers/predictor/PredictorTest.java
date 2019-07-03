/**
 *
 */
package org.theseed.genomes.kmers.predictor;

import java.io.IOException;
import org.theseed.genomes.kmers.DnaKmer;
import org.theseed.locations.Frame;
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

}
