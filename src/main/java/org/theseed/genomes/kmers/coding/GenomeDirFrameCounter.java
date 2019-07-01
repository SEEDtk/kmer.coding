/**
 *
 */
package org.theseed.genomes.kmers.coding;

import java.io.File;
import java.io.IOException;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genomes.GenomeDirectory;
import org.theseed.genomes.kmers.DnaKmer;

/**
 *
 * This is the primary class for computing the KmerFrameCounter from a directory of genomes.
 * The command-line options are
 *
 * 	K		kmer size (default is 15)
 *
 * The positional parameters are the name of the input directory (which must exist) and the
 * name of the output directory (which may need to be created).
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
    @Option(name="-K", aliases= {"--kmer"}, metaVar="15", usage="kmer size")
    private void setKmer(int newSize) {
        DnaKmer.setSize(newSize);
    }

    /** input directory name */
    @Argument(index=0, metaVar="inputDir", usage="input GTO directory",
            required=true, multiValued=false)
    private File inputDir;

    /** output directory */
    @Argument(index=1, metaVar="outDir", usage="output result directory",
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
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                this.inputGenomes = new GenomeDirectory(this.inputDir.getPath());
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
        System.err.println("Input directory is " + this.inputGenomes + ".");
        System.err.println("Output directory is " + this.outDir + ".");
        System.err.println("Kmer size is " + DnaKmer.getSize() + ".");
        // TODO Create the kmer counter.
        // TODO loop through the genomes, counting kmers
        // TODO output significant kmers
    }

}
