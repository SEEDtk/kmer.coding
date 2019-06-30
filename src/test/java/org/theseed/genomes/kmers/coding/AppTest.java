package org.theseed.genomes.kmers.coding;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.io.FileNotFoundException;

import com.github.cliftonlabs.json_simple.JsonException;

import org.theseed.genomes.Feature;
import org.theseed.genomes.Genome;
import org.theseed.genomes.Location;
import org.theseed.genomes.Frame;
import org.theseed.genomes.kmers.DnaKmer;
import org.theseed.genomes.kmers.SequenceDnaKmers;

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
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    private static final String myProtein = "MNERYQCLKTKEYQALLSSKGRQIFAKRKIDMKSVFGQIKVCLGYKRCHLRGKRQVRIDMGFILMANNLLKYNKRKRQN";
    private static final String myDna1 = "atgaatgaacgttaccagtgtttaaaaactaaagaatatcaggcacttttatcttccaagggtagacaaattttcgctaaacgtaagattgatatgaaatctgtctttgggcagataaaggtttgtttgggttataagagatgtcatctgagaggtaagcgtcaagtgagaattgacatgggattcatactcatggccaacaacctgctgaaatataataagagaaagaggcaaaattaa";
    private static final String myDna2 = "aaatagatttcaaaatgataaaaacgcatcctatcaggtttgagtgaacttgataggatgcgttttagaatgtcaaaattaattgagtttg";

    /**
     * Main test of genomes.
     * @throws FileNotFoundException
     * @throws JsonException
     */
    public void testGenome() throws FileNotFoundException, JsonException
    {

        Genome myGto = new Genome("gto_test/1313.7001.gto");
        assertEquals("Genome ID not correct.", "1313.7001", myGto.getId());
        assertEquals("Genome name not correct.", "Streptococcus pneumoniae P210774-233", myGto.getName());
        assertNull("Nonexistent feature found.", myGto.getFeature("fig|1313.7001.cds.75"));
        // Now we need to pull out a PEG and ask about it.
        Feature myFeature = myGto.getFeature("fig|1313.7001.peg.758");
        assertNotNull("Sample feature not found.", myFeature);
        assertEquals("Incorrect feature found.", "fig|1313.7001.peg.758", myFeature.getId());
        assertEquals("Incorrect function in sample feature.", "Transposase, IS4 family", myFeature.getFunction());
        assertEquals("Incorrect protein for sample feature.", myProtein, myFeature.getProteinTranslation());
        assertEquals("Incorrect DNA for sample feature.", myDna1, myGto.getDna("fig|1313.7001.peg.758"));
        // Next the location.
        Location myLoc = myFeature.getLocation();
        assertEquals("Incorrect contig for feature.", "1313.7001.con.0017", myLoc.getContigId());
        assertEquals("Incorrect left for feature.", 30698, myLoc.getLeft());
        assertEquals("Incorrect right for feature", 30937, myLoc.getRight());
        assertEquals("Incorrect begin for feature.", 30937, myLoc.getBegin());
        assertEquals("Incorrect length for feature.", 240, myLoc.getLength());
        assertEquals("Incorrect strand for feature.", '-', myLoc.getDir());
        assertFalse("Segmentation flag failure.", myLoc.isSegmented());
        // Now we check a segmented location.
        myFeature = myGto.getFeature("fig|1313.7001.repeat_unit.238");
        assertEquals("Incorrect DNA for segmented feature.", myDna2, myGto.getDna(myFeature.getId()));
        myLoc = myFeature.getLocation();
        assertEquals("Incorrect contig for segmented location.", "1313.7001.con.0018", myLoc.getContigId());
        assertEquals("Incorrect left for segmented location.", 11908, myLoc.getLeft());
        assertEquals("Incorrect right for segmented location.", 12116, myLoc.getRight());
        assertEquals("Incorrect begin for feature.", 11908, myLoc.getBegin());
        assertEquals("Incorrect length for feature.", 209, myLoc.getLength());
        assertEquals("Incorrect strand for segmented location.", '+', myLoc.getDir());
        assertTrue("Segmentation flag failure.", myLoc.isSegmented());
        // Now iterate over the proteins.
        for (Feature feat : myGto.getPegs()) {
            assertEquals("Feature" + feat.getId() + " is not a PEG.", "CDS", feat.getType());
        }
    }

    private static final String mySequence = "atgaatgaacgttaccagtgtttaaaaactxaaagaatatcaggcacttttatcttccaa";
                                           // 000000000111111111122222222223333333333444444444455555555555
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
        SequenceDnaKmers iterator = new SequenceDnaKmers(mySequence);
        while (iterator.nextKmer()) {
            int idx = iterator.idx();
            assertTrue("Invalid kmer returned at position " + iterator.getPos() + ".", idx >= 0);
            int pos = iterator.getPos();
            String kmer = iterator.toString();
            assertEquals("Iterated kmer did not recurse at position " + pos + ".",
                    mySequence.substring(pos-1, pos+9), kmer);
        }
        assertEquals("Iterator ended too soon.", mySequence.length() - 9, iterator.getPos());
    }

    /**
     * Main test of frames.
     */
    public void testFrames() {
        // Create a location for our main sequence.
        Location myLoc = Location.create("mySequence", "+");
        myLoc.addRegion(10, 20);
        myLoc.addRegion(30, 20);
        // Insure the kmer size is 15.
        DnaKmer.setSize(15);
        // Here are our expected results. Note there is no result for position 0.  We will start at 1.
        Frame[] results = new Frame[]
            { Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.P3, Frame.P1, Frame.P2, Frame.P3, Frame.P1, Frame.P2, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.P3, Frame.P1, Frame.P2, Frame.P3, Frame.P1, Frame.P2, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
              Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0
            };
        for (int pos = 1; pos < 60; pos++) {
            Frame computed = myLoc.kmerFrame(pos, DnaKmer.getSize());
            assertEquals("Incorrect frame result (plus strand " + pos + ").", results[pos], computed);
        }
        myLoc = Location.create("mySequence", "-");
        myLoc.addRegion(29, 20);
        myLoc.addRegion(49, 20);
        results = new Frame[]
                { Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.M2, Frame.M1, Frame.M3, Frame.M2, Frame.M1, Frame.M3, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.M2, Frame.M1, Frame.M3, Frame.M2, Frame.M1, Frame.M3, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX, Frame.XX,
                  Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0, Frame.F0
                };
        for (int pos = 1; pos < 60; pos++) {
            Frame computed = myLoc.kmerFrame(pos, DnaKmer.getSize());
            assertEquals("Incorrect frame result (minus strand " + pos + ").", results[pos], computed);
        }

    }

    /**
     * Main test of locations.
     */
    public void testLocations() {

    }
}

