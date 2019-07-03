/**
 *
 */
package org.theseed.genomes.kmers.coding;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collection;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genomes.Contig;
import org.theseed.genomes.Genome;
import org.theseed.genomes.GenomeDirectory;
import org.theseed.genomes.kmers.DnaKmer;
import org.theseed.genomes.kmers.SequenceDnaKmers;
import org.theseed.genomes.kmers.SequenceDnaNormalKmers;
import org.theseed.genomes.kmers.SequenceDnaSpacedKmers;
import org.theseed.genomes.kmers.predictor.FramePredictor;
import org.theseed.locations.Frame;
import org.theseed.locations.LocationList;
import org.theseed.utils.CountMap;

import com.github.cliftonlabs.json_simple.JsonException;

/**
 *
 * This is the primary class for computing the KmerFrameCounter from a directory of genomes.
 * The command-line options are
 *
 * 	-K		kmer size and type (default is 15); use a number for normal kmers, a number followed
 * 			by "p" for spaced kmers
 * 	-t		minimum best-fraction for a useful kmer (default is 0.80)
 * 	-m		minimum best-hits for a useful kmer (default is 30)
 * 	-i		input directory containing the genomes-- if omitted, a previously-built database is
 * 			loaded from the output directory
 *
 * The positional parameter is the name of the output directory (which may need to be created).
 *
 * @author Bruce Parrello
 */
public class GenomeDirFrameCounter {

    // FIELDS

    /** object to manage input directory */
    private GenomeDirectory inputGenomes;

    /** type of kmer processing */
    Class<? extends SequenceDnaKmers> kmerType;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** kmer size and type to use (default 15); suffix "p" indicates spaced kmers
     * @throws CmdLineException */
    @Option(name="-K", aliases={"--kmer"}, metaVar="15", usage="kmer size (XX for normal, XXp for spaced)")
    private void setKmer(String newSize) {
        int realSize;
        // Determine the kmer type.
        if (newSize.endsWith("p")) {
            kmerType = SequenceDnaSpacedKmers.class;
            realSize = Integer.valueOf(newSize.substring(0, newSize.length() - 1));
            // Insure the size is valid for a spaced kmer.
            if (realSize % 2 != 0) {
                throw new IllegalArgumentException("Spaced kmer sizes must be a multiple of 2.");
            }
        } else {
            kmerType = SequenceDnaNormalKmers.class;
            realSize = Integer.valueOf(newSize);
            // Insure the size is valid for a normal kmer.
            if (realSize % 3 != 0) {
                throw new IllegalArgumentException("Normal kmer sizes must be a multiple of 3.");
            }
        }
        // Store the kmer size.
        DnaKmer.setSize(realSize);
    }

    /** genome directory for optional testing set */
    @Option(name="-testDir", aliases={"--testDir"}, metaVar="genomeDir", usage="directory of GTOs for testing")
    private File testDir;

    /** minimum fraction for a kmer to be useful */
    @Option(name="-t", aliases={"--threshold"}, metaVar="0.8", usage="best-percent threshold for a useful kmer")
    private double threshold;

    /** minimum hits to the best frame for a kmer to be useful */
    @Option(name="-m", aliases= {"--minHits"}, metaVar="0", usage="minimum hits in the best frame for a useful kmer")
    private int minHits;

    /** input directory name; if omitted, the kmer database is reloaded from the output directory */
    @Option(name="-i", aliases={"--inputDir"}, metaVar="inputDir", usage="input GTO directory")
    private File inputDir;

    /** output directory */
    @Argument(index=0, metaVar="outDir", usage="output result directory",
            required=true, multiValued=false)
    private File outDir;

    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.threshold = 0.80;
        this.minHits = 30;
        this.inputDir = null;
        this.testDir = null;
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                if (this.inputDir != null) {
                    this.inputGenomes = new GenomeDirectory(this.inputDir.getPath());
                }
                if (this.outDir.isDirectory()) {
                    retVal = true;
                } else {
                    retVal = this.outDir.mkdirs();
                    if (! retVal) {
                        System.err.println("Error creating output directory" + this.outDir.getPath());
                    }
                }
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println("Error processing genome directory: " + e.getMessage());
        }
        return retVal;
    }

    public void run() {
        // Display the parameters.
        if (this.inputDir != null) {
            System.err.println("Input directory is " + this.inputGenomes + ".");
        }
        if (this.testDir != null) {
            System.err.println("Testing directory is " + this.testDir + ".");
        }
        System.err.println("Output directory is " + this.outDir + ".");
        System.err.println("Kmer size is " + DnaKmer.getSize() + ".");
        try {
            // Compute the kmer file name.
            File saveFile = new File(this.outDir, "kmers.ser");
            // Create the kmer counter.
            KmerFrameCounter bigCounter;
            if (this.inputDir != null) {
                // Here we have to create the kmer counter from the input directory.
                bigCounter = new KmerFrameCounter(this.kmerType);
                // Process the genomes.
                int gCount = 0;
                long start = System.currentTimeMillis();
                for (Genome genome : this.inputGenomes) {
                    gCount++;
                    System.err.println("Processing #" + gCount + ": " + genome + ".");
                    bigCounter.processGenome(genome);
                    // Display a time estimate every 100 genomes.
                    if (gCount % 100 == 0) {
                        double secsPerGenome = ((double) (System.currentTimeMillis() - start)) / (1000 * gCount);
                        double remainingMinutes = (this.inputGenomes.size() - gCount) * secsPerGenome / 60;
                        System.err.printf("TIME ESTIMATE: %4.2f seconds/genome, %4.1f minutes left.\n",
                                secsPerGenome, remainingMinutes);
                    }
                }
                System.err.println("Saving results.");
                bigCounter.save(saveFile);
            } else {
                // Here we have to reload an existing kmer counter database.
                System.err.println("Loading saved kmer database.");
                long start = System.currentTimeMillis();
                bigCounter = new KmerFrameCounter(saveFile);
                this.kmerType = bigCounter.getKmerType();
                double timeToLoad = ((double) (System.currentTimeMillis() - start)) / 1000;
                System.err.printf("%4.2f seconds to load database.\n", timeToLoad);
            }
            System.err.println("Searching for useful kmers.");
            // Open the kmer output file.
            File kmerFile = new File(this.outDir, "kmers.tbl");
            PrintWriter kmerWriter = new PrintWriter(kmerFile);
            // Start with a header.
            kmerWriter.println("kmer\tframe\tfraction\thits");
            // This will count the good kmers found.
            int goodCount = 0;
            // These are used to compute the mean hits and fraction.
            double totalFrac = 0.0;
            int countKmers = 0;
            long totalHits = 0;
            // This will count the good kmers per frame.
            int[] found = new int[Frame.nFrames];
            Arrays.fill(found, 0);
            long start = System.currentTimeMillis();
            // Loop through all the kmers.
            for (DnaKmer kmer : bigCounter) {
                Frame bestFrame = bigCounter.getBest(kmer);
                double frac = bigCounter.getFrac(kmer, bestFrame);
                int hits = bigCounter.getCount(kmer, bestFrame);
                totalFrac += frac;
                totalHits += hits;
                countKmers++;
                if (frac > this.threshold && hits > this.minHits) {
                    // Here the kmer is good enough.
                    kmerWriter.format("%s\t%s\t%04.2f\t%d%n", kmer, bestFrame, frac, hits);
                    goodCount++;
                    found[bestFrame.ordinal()]++;
                }
            }
            double secsToSearch = ((double) (System.currentTimeMillis() - start)) / 1000;
            System.err.format("%4.2f seconds to search kmer database%n", secsToSearch);
            kmerWriter.close();
            System.err.println("Writing report.");
            // Open the report output file.
            File reportFile = new File(this.outDir, "kmers.report.txt");
            PrintWriter reportWriter = new PrintWriter(reportFile);
            // Write the report.
            reportWriter.println(goodCount + " good kmers found.");
            reportWriter.println(countKmers + " unique kmers found.");
            double meanFrac = totalFrac / countKmers;
            double meanHits = ((double) totalHits) / countKmers;
            reportWriter.format("Mean hits per kmer = %4.2f, mean fraction = %4.2f%n",
                    meanHits, meanFrac);
            this.testFramePredictions(kmerFile, reportWriter, found);
            reportWriter.close();
            System.err.println("All done.");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Produce a report on the predictive power of a kmer set produced by this object.
     *
     * @param predFile		"kmers.tbl" file output by this object.
     * @param reportWriter	output writer for the report
     * @param found			an array of kmers found for each frame
     *
     * @throws NumberFormatException
     * @throws IOException
     * @throws JsonException
     */
    private void testFramePredictions(File predFile, PrintWriter reportWriter, int[] found)
            throws NumberFormatException, IOException, JsonException {
        // We will track counts in here.
        CountMap<DnaKmer> testCounts = new CountMap<DnaKmer>();
        int misses = 0;
        // Load the predictor.
        long start = System.currentTimeMillis();
        FramePredictor testPred = new FramePredictor(predFile.getPath());
        double loadTime = ((double) (System.currentTimeMillis() - start) / 1000);
        System.err.format("Predictor load test successful. %4.2f seconds.%n", loadTime);
        // Do we have genomes to test?
        if (this.testDir == null) {
            // No.  Write out the frame counts.
            reportWriter.format("%-8s %8s%n", "Frame", "kmers");
            for (Frame frm : Frame.all) {
                reportWriter.format("%-8s %8d %8d%n", frm, found[frm.ordinal()]);
            }
        } else {
            // Yes.  Load the genomes so we can test the predictor.
            GenomeDirectory genomes = new GenomeDirectory(this.testDir.getPath());
            for (Genome myGto : genomes) {
                System.err.println("Testing against " + myGto);
                Map<String, LocationList> gtoMap = LocationList.createGenomeCodingMap(myGto);
                // Loop through the contigs.
                Collection<Contig> allContigs = myGto.getContigs();
                for (Contig contig : allContigs) {
                    SequenceDnaKmers contigKmers = SequenceDnaKmers.build(this.kmerType, contig.getSequence());
                    LocationList contigLocs = gtoMap.get(contig.getId());
                    while (contigKmers.nextKmer()) {
                        int pos = contigKmers.getPos();
                        Frame predicted = testPred.frameOf(contigKmers);
                        if (predicted == Frame.XX) {
                            misses++;
                        } else {
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
            // Output what we found.
            double hitPercent = ((double) (goodHits + badHits) * 100) / (goodHits + badHits + misses);
            reportWriter.format("%d genomes were examined.%n", genomes.size());
            reportWriter.format("Total hits = %d good, %d bad. Percent hits %4.2f. Total misses = %d.%n",
                    goodHits, badHits, hitPercent, misses);
            reportWriter.format("%-8s %8s %8s %8s %8s%n", "Frame", "kmers", "goodHits", "badHits", "%good");
            for (Frame frm : Frame.sorted) {
                int good = frameCounts.good(frm);
                int bad = frameCounts.bad(frm);
                double goodPercent = (good <= 0 ? 0 : ((double) (good * 100)) / (good + bad));
                reportWriter.format("%-8s %8d %8d %8d %8.2f %n", frm, found[frm.ordinal()],
                        good, bad, goodPercent);
            }
        }
    }

    /**
     * @return the type of kmers used for counting
     */
    public Class<? extends SequenceDnaKmers> getKmerType() {
        return this.kmerType;
    }

    /**
     * @return the kmer size used for counting
     */
    public int getKmerSize() {
        return DnaKmer.getSize();
    }

    /**
     * @return the number of input genomes, or 0 if we are restoring from a saved file
     */
    public int getInputGenomesCount() {
        int retVal = 0;
        if (this.inputDir != null) {
            retVal = this.inputGenomes.size();
        }
        return retVal;
    }

    /**
     * @return the name of the test directory, or an empty string if there is none
     */
    public String getTestDir() {
        String retVal = "";
        if (this.testDir != null) {
            retVal = this.testDir.getPath();
        }
        return retVal;
    }


}