/**
 *
 */
package org.theseed.genome.kmers.predictor;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Scanner;

import org.theseed.genome.kmers.DnaKmer;
import org.theseed.locations.Frame;

/**
 * This class manages a data structure that maps DnaKmer objects to coding frames.
 *
 * @author Bruce Parrello
 *
 */
public class FramePredictor {

    // FIELDS

    private HashMap<DnaKmer, Frame> kmerMap;

    /**
     * Load a frame predictor from a file.  This file is "kmers.tbl" output by GenomeDirFrameCounter.
     *
     * @param fileName	name of the file containing the input table
     * @throws FileNotFoundException
     */
    public FramePredictor(String fileName) throws FileNotFoundException {
        // Create the kmer map.
        this.kmerMap = new HashMap<DnaKmer, Frame>();
        // Open the input file as a scanner.
        File inFile = new File(fileName);
        Scanner fileReader = new Scanner(inFile);
        // Throw away the header line.
        fileReader.nextLine();
        // Loop through the data lines.
        while (fileReader.hasNext()) {
            // Read the kmer and the target frame.
            DnaKmer kmer = new DnaKmer(fileReader.next());
            Frame bestFrame = Frame.frameOf(fileReader.next());
            this.kmerMap.put(kmer, bestFrame);
            // Skip the statistical data for now.  When we are more sophisticated, we will
            // use them to compute weights.
            fileReader.nextLine();
        }
        fileReader.close();
    }

    /**
     * @return the frame predicted by the specified kmer
     *
     * @param kmer kmer whose prediction is desired
     */
    public Frame frameOf(DnaKmer kmer) {
        return this.kmerMap.getOrDefault(kmer, Frame.XX);
    }

    /**
     * @return the frame predicted by the specified kmer
     *
     * @param kmer string form of kmer whose prediction is desired
     */
    public Frame frameOf(String kmerString) {
        DnaKmer kmer = new DnaKmer(kmerString);
        return this.frameOf(kmer);
    }


}
