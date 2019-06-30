package org.theseed.genomes.kmers.coding;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.io.FileNotFoundException;
import java.util.Collection;
import java.util.Iterator;

import com.github.cliftonlabs.json_simple.JsonException;

import org.theseed.genomes.Feature;
import org.theseed.genomes.Genome;
import org.theseed.genomes.Location;
import org.theseed.genomes.LocationList;
import org.theseed.genomes.Region;
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
     * Test of plus-locations.
     */
    public void testPlusLocations() {
        Location plusLoc = Location.create("mySequence1", "+");
        plusLoc.addRegion(402, 500);
        assertFalse("Segmentation error on single segment.", plusLoc.isSegmented());
        assertTrue("Location does not default to valid.", plusLoc.isValid());
        assertEquals("Invalid begin in plusLoc.", 402, plusLoc.getBegin());
        plusLoc.addRegion(10, 200);
        assertEquals("Invalid left position in plusLoc.", 10, plusLoc.getLeft());
        assertEquals("Invalid right position in plusLoc.", 901, plusLoc.getRight());
        assertEquals("Invalid contig in plusLoc.", "mySequence1", plusLoc.getContigId());
        assertEquals("Invalid strand in plusLoc.", '+', plusLoc.getDir());
        assertEquals("Invalid length in plusLoc.", 892, plusLoc.getLength());
        assertEquals("Invalid begin in plusLoc after add.", 10, plusLoc.getBegin());
        plusLoc.putRegion(1, 6);
        assertEquals("Invalid left position in plusLoc after put.", 1, plusLoc.getLeft());
        assertEquals("Invalid right position in plusLoc after put.", 901, plusLoc.getRight());
        plusLoc.putRegion(250, 349);
        assertEquals("Invalid left position in plusLoc after internal put.", 1, plusLoc.getLeft());
        assertEquals("Invalid right position in plusLoc after internal put.", 901, plusLoc.getRight());
        assertTrue("Segmentation error on multiple segments.", plusLoc.isSegmented());
        plusLoc.invalidate();
        assertFalse("Location did not invalidate.", plusLoc.isValid());
        plusLoc.setLeft(260);
        Collection<Region> plusRegions = plusLoc.getRegions();
        for (Region region : plusRegions) {
            int left = region.getLeft();
            assertTrue("Invalid region " + region + " extends past 260.", left >= 260);
        }
        Location cloneLoc = (Location) plusLoc.clone();
        assertEquals("Clone contig does not match.", plusLoc.getContigId(), cloneLoc.getContigId());
        assertEquals("Clone direction does not match.", plusLoc.getDir(), cloneLoc.getDir());
        assertEquals("Clone left does not match.", plusLoc.getLeft(), cloneLoc.getLeft());
        assertEquals("Clone length does not match.", plusLoc.getLength(), cloneLoc.getLength());
        assertTrue("Segmentation error on clone.", cloneLoc.isSegmented());
        assertFalse("Validation error on clone.", cloneLoc.isValid());
        Collection<Region> cloneRegions = cloneLoc.getRegions();
        for (Region region : plusRegions) {
            assertTrue("Region " + region + " not found in clone.", region.containedIn(cloneRegions));
        }
    }

    /**
     * Test of minus-locations.
     */
    public void testMinusLocations() {
        Location minusLoc = Location.create("mySequence1", "-");
        minusLoc.addRegion(901, 500);
        assertFalse("Segmentation error on single segment.", minusLoc.isSegmented());
        assertTrue("Location does not default to valid.", minusLoc.isValid());
        assertEquals("Invalid begin in minusLoc.", 901, minusLoc.getBegin());
        minusLoc.addRegion(209, 200);
        assertEquals("Invalid left position in minusLoc.", 10, minusLoc.getLeft());
        assertEquals("Invalid right position in minusLoc.", 901, minusLoc.getRight());
        assertEquals("Invalid contig in minusLoc.", "mySequence1", minusLoc.getContigId());
        assertEquals("Invalid strand in minusLoc.", '-', minusLoc.getDir());
        assertEquals("Invalid length in minusLoc.", 892, minusLoc.getLength());
        assertEquals("Invalid begin in minusLoc after add.", 901, minusLoc.getBegin());
        minusLoc.putRegion(1, 6);
        assertEquals("Invalid left position in minusLoc after put.", 1, minusLoc.getLeft());
        assertEquals("Invalid right position in minusLoc after put.", 901, minusLoc.getRight());
        minusLoc.putRegion(250, 349);
        assertEquals("Invalid left position in minusLoc after internal put.", 1, minusLoc.getLeft());
        assertEquals("Invalid right position in minusLoc after internal put.", 901, minusLoc.getRight());
        assertTrue("Segmentation error on multiple segments.", minusLoc.isSegmented());
        minusLoc.invalidate();
        assertFalse("Location did not invalidate.", minusLoc.isValid());
        minusLoc.setRight(260);
        Collection<Region> minusRegions = minusLoc.getRegions();
        for (Region region : minusRegions) {
            int right = region.getRight();
            assertTrue("Invalid region extends to " + right + " past 260.", right <= 260);
        }
        Location cloneLoc = (Location) minusLoc.clone();
        assertEquals("Clone contig does not match.", minusLoc.getContigId(), cloneLoc.getContigId());
        assertEquals("Clone direction does not match.", minusLoc.getDir(), cloneLoc.getDir());
        assertEquals("Clone left does not match.", minusLoc.getLeft(), cloneLoc.getLeft());
        assertEquals("Clone length does not match.", minusLoc.getLength(), cloneLoc.getLength());
        assertTrue("Segmentation error on clone.", cloneLoc.isSegmented());
        assertFalse("Validation error on clone.", cloneLoc.isValid());
        Collection<Region> cloneRegions = cloneLoc.getRegions();
        for (Region region :minusRegions) {
            assertTrue("Region " + region + " not found in clone.", region.containedIn(cloneRegions));
        }
    }

    /**
     * Location comparison test
     */
    public void testLocations() {
        Location loc1 = Location.create("myContig", "+", 1000, 1999);
        Location loc2 = Location.create("myContig", "-", 1100, 1199);
        Location loc3 = Location.create("myContig", "-", 1150, 1249);
        Location loc4 = Location.create("yourContig", "-", 1150, 1249);
        assertTrue("loc1 is not less than loc2.", loc1.compareTo(loc2) < 0);
        assertTrue("loc2 is not less than loc3.", loc2.compareTo(loc3) < 0);
        assertTrue("loc1 does not contain loc2.", loc1.contains(loc2));
        assertFalse("loc2 contains loc3.", loc2.contains(loc3));
        assertFalse("loc3 contains loc2.", loc3.contains(loc2));
        assertFalse("loc1 contains loc4.", loc1.contains(loc4));
    }

    /**
     * Main location list test
     */
    public void testLocationList() {
        LocationList newList = new LocationList("myContig");
        Location[] locs = { Location.create("myContig", "+", 10, 99, 102, 199),
                            Location.create("myContig", "-", 100, 400),
                            Location.create("myContig", "-", 500, 999),
                            Location.create("myContig", "-", 3000, 3199),
                            Location.create("myContig", "+", 4000, 4999),
                            Location.create("myContig", "-", 4100, 4199),
                            Location.create("myContig", "-", 6000, 6199),
                            Location.create("myContig", "+", 5000, 5099),
                            Location.create("myContig", "-", 5000, 5099),
                            Location.create("myContig", "+", 5200, 5299),
                            Location.create("myContig", "-", 5250, 5299),
                            Location.create("myContig", "+", 5350, 5399),
                            Location.create("myContig", "-", 5300, 5399),
                            Location.create("myContig", "+", 6800, 7299),
                            Location.create("myContig", "+", 8000, 8199),
                            Location.create("myContig", "+", 8400, 8599),
                            Location.create("myContig", "-", 8100, 8449),
                            Location.create("myContig", "+", 9100, 9199),
                            Location.create("myContig", "+", 9200, 9299),
                            Location.create("myContig", "-", 9300, 9399),
                            Location.create("myContig", "+", 9400, 9499),
                            Location.create("myContig", "-", 9000, 9999),
                          };
        // Add all these locations to the list.
        for (Location loc : locs) {
            newList.addLocation(loc);
        }
        // Now we need to verify that none of the stored locations overlap.
        Iterator<Location> iter = newList.iterator();
        Location prev = iter.next();
        while (iter.hasNext()) {
            Location next = iter.next();
            assertTrue(prev + " overlaps " + next, prev.getRight() < next.getLeft());
        }
        // Next, we want to show that every location in the location list is wholly
        // contained in one or more locations of the input list.  If it is contained
        // in more than one, or it is contained in a segmented one, it should be
        // invalid.  Otherwise, it should be on the same strand.
        for (Location loc : newList) {
            // This will remember the strand of the last location found.
            char strand = '0';
            // This will be set to TRUE if a location is segmented.
            Boolean invalid = false;
            // This will count the locations found.
            int found = 0;
            for (Location loc2 : locs) {
                if (loc2.contains(loc)) {
                    strand = loc2.getDir();
                    if (loc2.isSegmented()) {
                        invalid = true;
                    }
                    found++;
                }
            }
            if (found > 1) {
                assertFalse(loc + " is valid but it is in " + found + " input locations.", loc.isValid());
            } else if (invalid) {
                assertFalse(loc + " is valid but it is in a segmented input location.", loc.isValid());
            } else if (found == 0) {
                fail(loc + " was not found in the input locations.");
            } else {
                assertEquals(loc + " is not on the correct strand.", strand, loc.getDir());
            }
        }
        // Now for each input location, we need to show that every position is in a location
        // of the location list.
        for (Location loc : locs) {
            for (int pos = loc.getLeft(); pos <= loc.getRight(); pos++) {
                assertTrue("Position " + pos + " of location " + loc + " is not in the list.",
                        newList.computeStrand(pos) != '0');
            }
        }
    }
}
