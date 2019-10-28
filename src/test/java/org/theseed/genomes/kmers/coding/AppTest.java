package org.theseed.genomes.kmers.coding;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.genome.kmers.DnaKmer;
import org.theseed.genome.kmers.SequenceDnaKmers;
import org.theseed.genome.kmers.SequenceDnaNormalKmers;
import org.theseed.genome.kmers.SequenceDnaSpacedKmers;
import org.theseed.genome.kmers.coding.GenomeDirFrameCounter;
import org.theseed.genome.kmers.coding.KmerFrameCounter;
import org.theseed.genome.kmers.predictor.FramePredictor;
import org.theseed.locations.Frame;
import org.theseed.locations.Location;
import org.theseed.locations.Region;

/**
 * Unit test for simple App.
 */
public class AppTest
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     * @throws IOException
     * @throws NumberFormatException
     */
    public AppTest( String testName ) throws NumberFormatException, IOException
    {
        super( testName );
        this.myGto = new Genome(new File("src/test/gto_test", "/1313.7001.gto"));
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    private Genome myGto = null;


    private static final String mySequence = "atgaatgaacgttaccagtgtttaaaaactxaaagaatatcaggcacttttatcttccaa";
                                           // 000000000111111111122222222223333333333444444444455555555556
                                           // 123456789012345678901234567890123456789012345678901234567890

    /**
     * Main test of kmers.
     */
    public void testKmers() {
        // Insure the kmer size is 15.
        DnaKmer.setSize(15);
        DnaKmer kmer1 = new DnaKmer(mySequence, 10);
        assertTrue("Kmer index should be positive.", kmer1.idx() >= 0);
        assertEquals("Kmer index did not recurse.", mySequence.substring(9, 24), kmer1.toString());
        assertEquals("Kmer revcmp is incorrect.", "taaacactggtaacg", kmer1.toRString());
        kmer1 = new DnaKmer(mySequence, 20);
        assertEquals("Invalid handling of ambiguity character", DnaKmer.NULL, kmer1.idx());
        DnaKmer.setSize(10);
        kmer1 = new DnaKmer(mySequence, 20);
        assertTrue("Kmer length did not set properly.", kmer1.idx() >= 0);
        assertEquals("Kmer index did not recurse at new length.", mySequence.substring(19, 29), kmer1.toString());
        SequenceDnaKmers iterator = new SequenceDnaNormalKmers(mySequence);
        while (iterator.nextKmer()) {
            int idx = iterator.idx();
            assertTrue("Invalid kmer returned at position " + iterator.getPos() + ".", idx >= 0);
            int pos = iterator.getPos();
            String kmer = iterator.toString();
            assertEquals("Iterated kmer did not recurse at position " + pos + ".",
                    mySequence.substring(pos-1, pos+9), kmer);
        }
        assertEquals("Iterator ended too soon.", 52, iterator.getPos());
        iterator = new SequenceDnaSpacedKmers(mySequence);
        while (iterator.nextKmer()) {
            int idx = iterator.idx();
            assertTrue("Invalid kmer returned at position " + iterator.getPos() + ".", idx >= 0);
            int pos = iterator.getPos();
            String kmer = iterator.toString();
            for (int k = 0; k < 5; k++) {
                assertEquals("Error at position " + 2*k + " in spaced kmer for position " + pos + ".",
                        mySequence.substring(pos+3*k-1, pos+3*k), kmer.substring(2*k, 2*k+1));
            }
        }
        assertEquals("Iterator ended too soon.", 46, iterator.getPos());
        // Test kmer comparison.
        DnaKmer.setSize(15);
        DnaKmer kmerA = new DnaKmer("actccagcaagcatc");
        DnaKmer kmerB = new DnaKmer("actccagcaagcatc");
        DnaKmer kmerC = new DnaKmer("gatgcttgctggagt");
        assertTrue("Equal returns FALSE for same kmers.", kmerA.equals(kmerB));
        assertEquals("Hash codes not equal for same kmers.", kmerA.hashCode(), kmerB.hashCode());
        assertTrue("Equal not commutative for same kmers.", kmerB.equals(kmerA));
        assertFalse("Equal returns TRUE for revcmp kmers.", kmerA.equals(kmerC));
        assertFalse("Hash codes equal for same revcmp kmers.", kmerA.hashCode() == kmerC.hashCode());
        assertEquals("Compare returns nonzero for same kmers.", 0, kmerA.compareTo(kmerB));
        assertFalse("Compare returns zero for revcmp kmers.", kmerA.compareTo(kmerC) == 0);
        assertEquals("Comparison not an ordering.", kmerA.compareTo(kmerC), -kmerC.compareTo(kmerA));
        // Test fromRString
        kmerA.setIdx(DnaKmer.fromRString("actccagcaagcatc"));
        assertEquals("Rstring did not produce proper kmer.", kmerA, kmerC);
        assertEquals("Short kmer failed in RString.", DnaKmer.EOF, DnaKmer.fromRString("aaaac"));
        assertEquals("Bad char failed in RString.", DnaKmer.NULL, DnaKmer.fromRString("acgtaxgtacgtcac"));
    }

    /**
     * Main test of frames.
     */
    public void testFrames() {
        // Verify frame transformations.
        assertEquals("P0 did not reverse.", Frame.M0, Frame.P0.rev());
        assertEquals("P1 did not reverse.", Frame.M1, Frame.P1.rev());
        assertEquals("P2 did not reverse.", Frame.M2, Frame.P2.rev());
        assertEquals("M0 did not reverse.", Frame.P0, Frame.M0.rev());
        assertEquals("M1 did not reverse.", Frame.P1, Frame.M1.rev());
        assertEquals("M2 did not reverse.", Frame.P2, Frame.M2.rev());
        assertEquals("F0 did not reverse.", Frame.F0, Frame.F0.rev());
        // Verify frame construction from labels.
        for (Frame frm : Frame.all) {
            assertEquals("Wrong frame for label " + frm, frm, Frame.frameOf(frm.toString()));
        }
        // Create a location for our main sequence.
        Location myLoc = Location.create("mySequence", "+");
        myLoc.addRegion(10, 20);
        myLoc.addRegion(30, 20);
        // Insure the kmer size is 15.
        DnaKmer.setSize(15);
        // Here are our expected results. Note there is no result for position 0.  We will start at 1.
        Frame[] results = new Frame[]
            { Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.P0, Frame.P1, Frame.P2, Frame.P0, Frame.P1, Frame.P2, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.P0, Frame.P1, Frame.P2, Frame.P0, Frame.P1, Frame.P2, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0
            };
        for (int pos = 1; pos < 60; pos++) {
            Frame computed = myLoc.regionFrame(pos, pos + DnaKmer.getSize() - 1);
            assertEquals("Incorrect frame result (plus strand " + pos + ").", results[pos], computed);
        }
        myLoc = Location.create("mySequence", "-");
        myLoc.addRegion(29, 20);
        myLoc.addRegion(49, 20);
        results = new Frame[]
                { Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.M2, Frame.M1, Frame.M0, Frame.M2, Frame.M1, Frame.M0, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.M2, Frame.M1, Frame.M0, Frame.M2, Frame.M1, Frame.M0, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0
                };
        for (int pos = 1; pos < 60; pos++) {
            Frame computed = myLoc.regionFrame(pos, pos + DnaKmer.getSize() - 1);
            assertEquals("Incorrect frame result (minus strand " + pos + ").", results[pos], computed);
        }

    }

    /**
     * The test for the insanely big kmer counter.
     *
     * @throws FileNotFoundException
     */
    public void testKmerCounter() throws FileNotFoundException {
        // Use small kmers to fit in debugger.
        DnaKmer.setSize(9);
        KmerFrameCounter bigCounter = new KmerFrameCounter(SequenceDnaNormalKmers.class);
        // Choose a kmer to manipulate.
        DnaKmer myKmer = new DnaKmer("actgtccat");
        // Verify that an uncounted kmer's best frame is invalid.
        assertEquals("Valid frame returned for empty kmer.", Frame.XX, bigCounter.getBest(myKmer));
        assertEquals("Nonzero fraction returned for empty kmer.", 0.0,
                bigCounter.getFrac(myKmer, Frame.F0), 0.0);
        // Test iterating the counts.
        for (Frame frm : Frame.all) {
            bigCounter.increment(myKmer, frm);
        }
        // Put some counts in.
        for (int i = 1; i <= 39; i++) {
            bigCounter.increment(myKmer, Frame.M0);
        }
        for (int i = 1; i <= 59; i++) {
            bigCounter.increment(myKmer, Frame.M2);
        }
        for (int i = 1; i <= 99; i++) {
            bigCounter.increment(myKmer, Frame.M1);
        }
        for (int i = 1; i <= 199; i++) {
            bigCounter.increment(myKmer, Frame.F0);
        }
        for (int i = 1; i <= 499; i++) {
            bigCounter.increment(myKmer, Frame.P1);
        }
        for (int i = 1; i <= 19; i++) {
            bigCounter.increment(myKmer, Frame.P2);
        }
        for (int i = 1; i <= 79; i++) {
            bigCounter.increment(myKmer, Frame.P0);
        }
        // Verify that the invalid frame doesn't crash us.
        bigCounter.increment(myKmer, Frame.XX);
        assertEquals("Invalid frame has nonzero count.", 0, bigCounter.getCount(myKmer, Frame.XX));
        // Verify the above counts.
        assertEquals("Count error in frame M3.", 40, bigCounter.getCount(myKmer, Frame.M0));
        assertEquals("Count error in frame M2.", 60, bigCounter.getCount(myKmer, Frame.M2));
        assertEquals("Count error in frame M1.", 100, bigCounter.getCount(myKmer, Frame.M1));
        assertEquals("Count error in frame F0.", 200, bigCounter.getCount(myKmer, Frame.F0));
        assertEquals("Count error in frame P1.", 500, bigCounter.getCount(myKmer, Frame.P1));
        assertEquals("Count error in frame P2.", 20, bigCounter.getCount(myKmer, Frame.P2));
        assertEquals("Count error in frame P3.", 80, bigCounter.getCount(myKmer, Frame.P0));
        // Verify the fractions.
        assertEquals("Frac error in frame M3.", 0.040, bigCounter.getFrac(myKmer, Frame.M0));
        assertEquals("Frac error in frame M2.", 0.060, bigCounter.getFrac(myKmer, Frame.M2));
        assertEquals("Frac error in frame M1.", 0.100, bigCounter.getFrac(myKmer, Frame.M1));
        assertEquals("Frac error in frame F0.", 0.200, bigCounter.getFrac(myKmer, Frame.F0));
        assertEquals("Frac error in frame P1.", 0.500, bigCounter.getFrac(myKmer, Frame.P1));
        assertEquals("Frac error in frame P2.", 0.020, bigCounter.getFrac(myKmer, Frame.P2));
        assertEquals("Frac error in frame P3.", 0.080, bigCounter.getFrac(myKmer, Frame.P0));
        // Verify the best frame.
        Frame bestFrame = bigCounter.getBest(myKmer);
        assertSame("Incorrect best frame.", Frame.P1, bestFrame);
        double frac = bigCounter.getFrac(myKmer, bestFrame);
        int hits = bigCounter.getCount(myKmer, bestFrame);
        // Test printing the best frame.
        PrintWriter kmerWriter = new PrintWriter("src/test/testOut.txt");
        kmerWriter.printf("%s\t%s\t%04.2f\t%d\n", myKmer, bestFrame, frac, hits);
        kmerWriter.close();
        // Check overflow.
        for (int i = 1; i <= 40000; i++) {
            bigCounter.increment(myKmer, Frame.P1);
        }
        assertEquals("Overflow occurred above 32K.", 40500, bigCounter.getCount(myKmer, Frame.P1));
        assertEquals("Invalid fraction at overflow point.", 0.988, bigCounter.getFrac(myKmer, Frame.P1), 0.001);
        assertEquals("Invalid frame choice at overflow point.", Frame.P1, bigCounter.getBest(myKmer));
        // Verify that only our one kmer is found by the iterator.
        for (DnaKmer kmer : bigCounter) {
            assertEquals("Incorrect kmer found.", kmer, myKmer);
        }
        // Erase the counter and add in a genome.
        bigCounter.clear();
        bigCounter.processGenome(this.myGto);
        // Checking the genome is hard.  We are going to check for "tgaatgaac"
        // in frame M1. We know this is in the M1 frame for one protein.
        DnaKmer targetKmer = new DnaKmer("gttcattca");
        assertTrue("Target kmer not in known frame.", bigCounter.getCount(targetKmer, Frame.M1) > 0);
        // Now we need to very serialization.  Write out the kmer counter.
        bigCounter.save("src/test/kmerTest.ser");
        // Save some key values.  We can't save them all, due to memory issues.
        int saveMyP3 = bigCounter.getCount(myKmer, Frame.P0);
        int saveMyM1 = bigCounter.getCount(myKmer, Frame.M1);
        int saveTargetP2 = bigCounter.getCount(targetKmer, Frame.P2);
        int saveTargetF0 = bigCounter.getCount(targetKmer, Frame.F0);
        // Erase the old kmer counter.
        bigCounter = null;
        // Read in the saved one.
        bigCounter = new KmerFrameCounter("src/test/kmerTest.ser");
        // Test the saved values.
        assertEquals("Error in saveMyP3.", saveMyP3, bigCounter.getCount(myKmer, Frame.P0));
        assertEquals("Error in saveMyM1.", saveMyM1, bigCounter.getCount(myKmer, Frame.M1));
        assertEquals("Error in saveTargetP2.", saveTargetP2, bigCounter.getCount(targetKmer, Frame.P2));
        assertEquals("Error in saveTargetF0.", saveTargetF0, bigCounter.getCount(targetKmer, Frame.F0));
    }

    /**
     * Test counter on spaced kmers.
     */
    public void testSpacedCounting() {
        DnaKmer.setSize(10);
        // Do some sanity checks on spaced reads.
        SequenceDnaSpacedKmers newProcessor = new SequenceDnaSpacedKmers("atgaatgaacgttaccagtgtttaaaaactaaagaatatcaggcacttttatct");
        assertTrue("Failed on first kmer read.", newProcessor.nextKmer());
        assertEquals("First kmer is wrong.", "ataagacgta", newProcessor.toString());
        newProcessor.reverse();
        assertEquals("Reverse of first kmer is wrong.", "gtacttatca", newProcessor.toString());
        assertTrue("Failed on second kmer read.", newProcessor.nextKmer());
        assertEquals("Second kmer is wrong.", "tgataagtac", newProcessor.toString());
        // Try a spaced reading of a genome.
        KmerFrameCounter bigCounter = new KmerFrameCounter(SequenceDnaSpacedKmers.class);
        bigCounter.processGenome(myGto);
        // Get the DNA for our sample peg.
        String dna = myGto.getDna("fig|1313.7001.peg.758");
        assertEquals("Wrong DNA for peg somehow.", "atgaatgaacgttac", dna.substring(0, 15));
        // Check it against the contig.
        Contig contig = myGto.getContig("1313.7001.con.0017");
        String origDna = contig.getDna(new Region(30923, 30937));
        assertEquals("Wrong DNA in contig somehow.", "gtaacgttcattcat", origDna);
        // Insure we found stuff we know is in our sample peg.
        // atgaatgaacgttaccagtgtttaaaaactaaagaatatcaggcacttttatct
        // at aa ga cg ta
        DnaKmer kmer = new DnaKmer("ataagacgta");
        assertTrue("ataagacgta not found in P0.", bigCounter.getCount(kmer, Frame.P0) > 0);
        kmer = new DnaKmer("tgataagtac");
        assertTrue("tgataagtac not found in P1.", bigCounter.getCount(kmer, Frame.P1) > 0);
        kmer = new DnaKmer("gatgacttcc");
        assertTrue("gatgacttcc not found in P2.", bigCounter.getCount(kmer, Frame.P2) > 0);
    }

    /**
     * Test command line
     */
    public void testCommandLine() {
        String[] args1 = {"--inputDir", "src/test/gto_test", "-K",  "8p",  "--testDir",
                "Bins_HMP", "TestKmers" };
        GenomeDirFrameCounter runObject = new GenomeDirFrameCounter();
        runObject.parseCommand(args1);
        assertEquals("Incorrect kmer type (8p).", SequenceDnaSpacedKmers.class, runObject.getKmerType());
        assertEquals("Incorrect kmer size (8p).", 8, runObject.getKmerSize());
        assertEquals("Incorrect input directory.", 4, runObject.getInputGenomesCount());
        assertEquals("Incorrect test directory.", "Bins_HMP", runObject.getTestDir());
        String args2[] = { "-K", "12", "TestKmers2" };
        runObject.parseCommand(args2);
        assertEquals("Incorrect kmer type (12).", SequenceDnaNormalKmers.class, runObject.getKmerType());
        assertEquals("Incorrect kmer size (12).", 12, runObject.getKmerSize());
        assertEquals("Incorrect use of input directory.", 0, runObject.getInputGenomesCount());
        assertEquals("Incorrect use of test directory.", "", runObject.getTestDir());
    }

    /**
     * Tests for FramePredictor
     * @throws IOException
     */
    public void testFramePredictor() throws IOException {
        FramePredictor testPred = new FramePredictor("src/test/kmers.tbl");
        DnaKmer testKmer = new DnaKmer("gacgggcgtgtagac");
        assertEquals("Test kmer in wrong frame.", Frame.M0, testPred.frameOf(testKmer));
        assertEquals("Non-coding kmer in wrong frame.", Frame.F0, testPred.frameOf("gacgggcggtgtgtg"));
        assertEquals("Plus-one kmer in wrong frame.", Frame.P1, testPred.frameOf("gacgggctacacatt"));
    }

}


