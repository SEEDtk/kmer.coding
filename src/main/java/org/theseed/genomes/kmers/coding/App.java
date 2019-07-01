package org.theseed.genomes.kmers.coding;

/**
 *
 * Process GTOs in a directory and count their kmers.
 *
 * @author Bruce Parrello
 */
public class App
{


    public static void main( String[] args )
    {
        GenomeDirFrameCounter runObject = new GenomeDirFrameCounter();
        boolean ok = runObject.parseCommand(args);
        if (ok) {
            runObject.run();
        }

    }
}
