/**
 *
 */
package org.theseed.genomes.kmers.coding;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genomes.Genome;
import org.theseed.genomes.GenomeDirectory;
import org.theseed.genomes.kmers.DnaKmer;
import org.theseed.locations.Frame;

/**
 *
 * This is the primary class for computing the KmerFrameCounter from a directory of genomes.
 * The command-line options are
 *
 * 	K		kmer size (default is 15)
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

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** kmer size to use (default 15) */
    @Option(name="-K", aliases={"--kmer"}, metaVar="15", usage="kmer size")
    private void setKmer(int newSize) {
        DnaKmer.setSize(newSize);
    }

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
        System.err.println("Output directory is " + this.outDir + ".");
        System.err.println("Kmer size is " + DnaKmer.getSize() + ".");
        try {
            // Compute the kmer file name.
            File saveFile = new File(this.outDir, "kmers.ser");
            // Create the kmer counter.
            KmerFrameCounter bigCounter;
            if (this.inputDir != null) {
                // Here we have to create the kmer counter from the input directory.
                bigCounter = new KmerFrameCounter();
                // Process the genomes.
                int gCount = 0;
                long start = System.currentTimeMillis();
                for (Genome genome : this.inputGenomes) {
                    gCount++;
                    System.err.println("Processing #" + gCount + ": " + genome + ".");
                    bigCounter.processGenome(genome);
                    // Display a time estimate every 100 genomes.
                    if (gCount % 100 == 0) {
                        double secsPerGenome = (System.currentTimeMillis() - start) / (1000 * gCount);
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
                bigCounter = KmerFrameCounter.load(saveFile);
                double timeToLoad = (System.currentTimeMillis() - start) / 1000;
                System.err.printf("%4.2f seconds to load database.\n", timeToLoad);
            }
            System.err.println("Searching for useful kmers.");
            // Open the kmer output file.
            File kmerFile = new File(this.outDir, "kmers.tbl");
            PrintWriter kmerWriter = new PrintWriter(kmerFile);
            // Start with a header.
            kmerWriter.print("kmer\tframe\tfraction\thits\n");
            // This will count the good kmers found.
            int goodCount = 0;
            // These are used to compute the mean hits and fraction.
            double totalFrac = 0.0;
            int countKmers = 0;
            long totalHits = 0;
            // This will count the good kmers per frame.
            int[] found = new int[Frame.nFrames];
            Arrays.fill(found, 0);
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
                    kmerWriter.printf("%s\t%s\t%04.2f\t%d\n", kmer, bestFrame, frac, hits);
                    goodCount++;
                    found[bestFrame.ordinal()]++;
                }
            }
            kmerWriter.close();
            System.err.println(goodCount + " good kmers found.");
            System.err.println(countKmers + " unique kmers found.");
            double meanFrac = totalFrac / countKmers;
            double meanHits = ((double) totalHits) / countKmers;
            System.err.printf("Mean hits per kmer = %4.2f, mean fraction = %4.2f", meanHits, meanFrac);
            System.err.println("Useful kmers per frame:");
            for (Frame frm : Frame.all) {
                System.err.printf("%-3s  %8d\n", frm, found[frm.ordinal()]);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
