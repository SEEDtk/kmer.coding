package org.theseed.genomes.kmers.coding;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import com.github.cliftonlabs.json_simple.JsonException;

import org.theseed.genomes.Feature;
import org.theseed.genomes.Genome;
import org.theseed.genomes.GenomeDirectory;
import org.theseed.genomes.kmers.DnaKmer;
import org.theseed.genomes.kmers.SequenceDnaKmers;
import org.theseed.genomes.kmers.SequenceDnaNormalKmers;
import org.theseed.genomes.kmers.SequenceDnaSpacedKmers;
import org.theseed.locations.Frame;
import org.theseed.locations.Location;
import org.theseed.locations.LocationList;
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
     * @throws JsonException
     * @throws FileNotFoundException
     * @throws NumberFormatException
     */
    public AppTest( String testName ) throws NumberFormatException, FileNotFoundException, JsonException
    {
        super( testName );
        this.myGto = new Genome("gto_test/1313.7001.gto");
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    private Genome myGto = null;


    private static final String myProtein = "MNERYQCLKTKEYQALLSSKGRQIFAKRKIDMKSVFGQIKVCLGYKRCHLRGKRQVRIDMGFILMANNLLKYNKRKRQN";
    private static final String myDna1 = "atgaatgaacgttaccagtgtttaaaaactaaagaatatcaggcacttttatcttccaagggtagacaaattttcgctaaacgtaagattgatatgaaatctgtctttgggcagataaaggtttgtttgggttataagagatgtcatctgagaggtaagcgtcaagtgagaattgacatgggattcatactcatggccaacaacctgctgaaatataataagagaaagaggcaaaattaa";
    private static final String myDna2 = "aaatagatttcaaaatgataaaaacgcatcctatcaggtttgagtgaacttgataggatgcgttttagaatgtcaaaattaattgagtttg";


    /**
     * Main test of genomes.
     * @throws FileNotFoundException
     * @throws JsonException
     */
    public void testGenome()
    {
        assertEquals("Genome ID not correct.", "1313.7001", this.myGto.getId());
        assertEquals("Genome name not correct.", "Streptococcus pneumoniae P210774-233", this.myGto.getName());
        assertNull("Nonexistent feature found.", this.myGto.getFeature("fig|1313.7001.cds.75"));
        // Now we need to pull out a PEG and ask about it.
        Feature myFeature = this.myGto.getFeature("fig|1313.7001.peg.758");
        assertNotNull("Sample feature not found.", myFeature);
        assertEquals("Incorrect feature found.", "fig|1313.7001.peg.758", myFeature.getId());
        assertEquals("Incorrect function in sample feature.", "Transposase, IS4 family", myFeature.getFunction());
        assertEquals("Incorrect protein for sample feature.", myProtein, myFeature.getProteinTranslation());
        assertEquals("Incorrect DNA for sample feature.", myDna1, this.myGto.getDna("fig|1313.7001.peg.758"));
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
        myFeature = this.myGto.getFeature("fig|1313.7001.repeat_unit.238");
        assertEquals("Incorrect DNA for segmented feature.", myDna2, this.myGto.getDna(myFeature.getId()));
        myLoc = myFeature.getLocation();
        assertEquals("Incorrect contig for segmented location.", "1313.7001.con.0018", myLoc.getContigId());
        assertEquals("Incorrect left for segmented location.", 11908, myLoc.getLeft());
        assertEquals("Incorrect right for segmented location.", 12116, myLoc.getRight());
        assertEquals("Incorrect begin for feature.", 11908, myLoc.getBegin());
        assertEquals("Incorrect length for feature.", 209, myLoc.getLength());
        assertEquals("Incorrect strand for segmented location.", '+', myLoc.getDir());
        assertTrue("Segmentation flag failure.", myLoc.isSegmented());
        // Now iterate over the proteins.
        for (Feature feat : this.myGto.getPegs()) {
            assertEquals("Feature" + feat.getId() + " is not a PEG.", "CDS", feat.getType());
        }
    }

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
        assertEquals("Iterator ended too soon.", 47, iterator.getPos());
        // Test kmer comparison.
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
    }

    /**
     * Main test of frames.
     */
    public void testFrames() {
        // Verify frame transformations.
        assertEquals("P3 did not reverse.", Frame.M0, Frame.P0.rev());
        assertEquals("P2 did not reverse.", Frame.M1, Frame.P2.rev());
        assertEquals("P1 did not reverse.", Frame.M2, Frame.P1.rev());
        assertEquals("M3 did not reverse.", Frame.P0, Frame.M0.rev());
        assertEquals("M2 did not reverse.", Frame.P1, Frame.M2.rev());
        assertEquals("M1 did not reverse.", Frame.P2, Frame.M1.rev());
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
        Region region1 = new Region(1000, 2000);
        Region region2 = new Region(1000, 2000);
        Region region3 = new Region(2000, 3000);
        assertTrue("region1 not equal to region 2.", region1.equals(region2));
        assertFalse("region1 equal to region 3.", region1.equals(region3));
        assertTrue("Region equals is not commutative.", region2.equals(region1));
        assertEquals("Equal regions have different hash codes.", region1.hashCode(), region2.hashCode());
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
        Location loc5 = Location.create("myContig", "+",  1000, 2000, 3000, 4000);
        Location loc6 = Location.create("myContig", "+", 3000, 4000, 1000, 2000);
        Location loc7 = Location.create("yourContig", "+", 1000, 2000, 3000, 4000);
        Location loc8 = Location.create("myContig", "-", 1000, 2000, 3000, 4000);
        Location loc9 = Location.create("myContig",  "+",  1000, 1999, 3000, 4000);
        assertTrue("loc5 not equal to loc6.", loc5.equals(loc6));
        assertEquals("Equal locations have different hash codes.", loc5.hashCode(), loc6.hashCode());
        assertFalse("Different contigs compare equal.", loc5.equals(loc7));
        assertFalse("Different strands compare equal.", loc5.equals(loc8));
        assertFalse("Different region counts compare equal.", loc5.equals(loc1));
        assertFalse("Different region extents compare equal.", loc5.equals(loc9));
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
            assertTrue("Failed to add " + loc + " to list.", newList.addLocation(loc));
        }
        Location badLoc = Location.create("yourContig", "+", 1000, 4999);
        assertFalse("Added wrong contig successfully.", newList.addLocation(badLoc));
        assertEquals("Invalid contig ID in location list.", "myContig", newList.getContigId());
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
            boolean invalid = false;
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
        // Finally, we want to test the frame computation.
        assertEquals("Invalid frame for pre-loc position.", Frame.XX, newList.computeRegionFrame(1, 15));
        assertEquals("Invalid frame for segmented position.", Frame.XX, newList.computeRegionFrame(40, 45));
        assertEquals("Invalid frame for simple minus position.", Frame.M1, newList.computeRegionFrame(390, 399));
        assertEquals("Invalid frame for simple plus position.", Frame.P0, newList.computeRegionFrame(4009, 4054));
        assertEquals("Invalid frame for near-overlap position.", Frame.P2, newList.computeRegionFrame(5235, 5245));
        assertEquals("Invalid frame for overlap position.", Frame.XX, newList.computeRegionFrame(5235, 5255));
        assertEquals("Invalid frame for extron position.", Frame.F0, newList.computeRegionFrame(7306, 7316));
    }

    /**
     * Basic test for location-list maps.
     */
    public void testContigMapping() {
        Map<String, LocationList> gList = LocationList.createGenomeCodingMap(this.myGto);
        LocationList contig0036 = gList.get("1313.7001.con.0036");
        assertNotNull("Contig 0036 not found.", contig0036);
        assertEquals("Incorrect strand found for test position 33996.", '+', contig0036.computeStrand(33996));
        assertEquals("Incorrect strand found for test position 30980.", '-', contig0036.computeStrand(30980));
        assertEquals("Incorrect strand found for test position 30984.", '0', contig0036.computeStrand(30984));
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
        PrintWriter kmerWriter = new PrintWriter("testOut.txt");
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
        DnaKmer targetKmer = new DnaKmer("tgaatgaac");
        assertTrue("Target kmer not in known frame.", bigCounter.getCount(targetKmer, Frame.M1) > 0);
        // Now we need to very serialization.  Write out the kmer counter.
        bigCounter.save("kmerTest.ser");
        // Save some key values.  We can't save them all, due to memory issues.
        int saveMyP3 = bigCounter.getCount(myKmer, Frame.P0);
        int saveMyM1 = bigCounter.getCount(myKmer, Frame.M1);
        int saveTargetP2 = bigCounter.getCount(targetKmer, Frame.P2);
        int saveTargetF0 = bigCounter.getCount(targetKmer, Frame.F0);
        // Erase the old kmer counter.
        bigCounter = null;
        // Read in the saved one.
        bigCounter = new KmerFrameCounter("kmerTest.ser");
        // Test the saved values.
        assertEquals("Error in saveMyP3.", saveMyP3, bigCounter.getCount(myKmer, Frame.P0));
        assertEquals("Error in saveMyM1.", saveMyM1, bigCounter.getCount(myKmer, Frame.M1));
        assertEquals("Error in saveTargetP2.", saveTargetP2, bigCounter.getCount(targetKmer, Frame.P2));
        assertEquals("Error in saveTargetF0.", saveTargetF0, bigCounter.getCount(targetKmer, Frame.F0));
    }

    /**
     * Test genome directories
     *
     * @throws IOException
     */
    public void testGenomeDir() throws IOException {
        GenomeDirectory gDir = new GenomeDirectory("gto_test");
        assertEquals("Wrong number of genomes found.", 4, gDir.size());
        // Run through an iterator.  We know the genome IDs, we just need to find them in order.
        String[] expected = new String[] { "1005394.4", "1313.7001", "1313.7002", "1313.7016" };
        int i = 0;
        for (Genome genome : gDir) {
            assertEquals("Incorrect result for genome at position " + i + ".", expected[i], genome.getId());
            i++;
        }
    }

}


