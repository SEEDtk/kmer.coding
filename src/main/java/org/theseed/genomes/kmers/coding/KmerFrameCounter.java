/**
 *
 */
package org.theseed.genomes.kmers.coding;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;

import org.theseed.genomes.Contig;
import org.theseed.genomes.Genome;
import org.theseed.genomes.kmers.DnaKmer;
import org.theseed.genomes.kmers.SequenceDnaKmers;
import org.theseed.locations.Frame;
import org.theseed.locations.LocationList;

/**
 * @author Bruce Parrello
 *
 * This class manages the frame counts for all the kmers.  It is a set of giant arrays, one per
 * frame, each array indexed by DnaKmer numbers.  To save space, the arrays contain short
 * integers, but these are treated as unsigned via some fancy numeric dancing.  All the
 * counts are returned as full integers.
 *
 */
public class KmerFrameCounter implements Iterable<DnaKmer> {

    // FIELDS
    /** the master array, indexed by frame ordinal and then kmer index */
    private short[][] countArray;
    /** the number of kmer values */
    private int size;


    /**
     * Construct an empty kmer frame counter.
     */
    public KmerFrameCounter() {
        this.size = DnaKmer.maxKmers();
        this.countArray = new short[7][this.size];
        this.clear();
    }

    /**
     * @return an iterator for all the kmers with nonzero counts
     */
    @Override
    public Iterator<DnaKmer> iterator() {
        return new KmerCountIterator();
    }

    /**
     * @return the kmer count for a specified frame.
     *
     * @param kmer	the relevant kmer
     * @param frm	the frame whose count is desired
     */
    public int getCount(DnaKmer kmer, Frame frm) {
        return iCount(kmer, frm.ordinal());
    }

    /**
     * @return the count for the frame with the specified ordinal
     *
     * @param kmer		the relevant kmer
     * @param ordinal	of the frame whose count is desired
     */
    private int iCount(DnaKmer kmer, int ordinal) {
        int retVal = 0;
        if (ordinal < 7) {
            retVal = (this.countArray[ordinal][kmer.idx()]) & 0xFFFF;
        }
        return retVal;
    }

    /**
     * Increment the kmer count for a specified frame.
     *
     * @param kmer	the relevant kmer
     * @param frm	the frame whose count is to be increments
     */
    public void increment(DnaKmer kmer, Frame frm) {
        if (frm != Frame.XX) {
            this.countArray[frm.ordinal()][kmer.idx()]++;
        }
    }

    /**
     * @return 	the best frame for a kmer; that is, the frame with the highest count,
     * 			or Frame.XX if the kmer has no instances
     *
     * @param kmer	the relevant kmer
     */
    public Frame getBest(DnaKmer kmer) {
        int best = 7;
        int bestCount = 0;
        for (int i = 0; i < 7; i++) {
            int count = this.iCount(kmer, i);
            if (count > bestCount) {
                bestCount = count;
                best = i;
            }
        }
        return Frame.idxFrame(best);
    }

    /**
     * @return	the fraction of time that the specified frame contains a kmer
     *
     * @param kmer	the relevant kmer
     */
    public double getFrac(DnaKmer kmer, Frame frm) {
        double retVal = 0;
        int frmCount = this.getCount(kmer, frm);
        int total = 0;
        for (int i = 0; i < 7; i++) {
            total += this.iCount(kmer, i);
        }
        if (total > 0) {
            retVal = ((double) frmCount) / ((double) total);
        }
        return retVal;
    }


    /**
     * Iterator class for finding kmers with nonzero counts.
     */
    public class KmerCountIterator implements Iterator<DnaKmer> {

        // FIELDS
        /** This always points to the next kmer to return.  If it is too high, we are at the end. */
        DnaKmer nextKmer;
        /** To avoid storage madness, this is the kmer we return to the user. */
        DnaKmer thisKmer;

        public KmerCountIterator() {
            this.thisKmer = new DnaKmer();
            // Find the first kmer (if any).
            this.findAfter(-1);
        }

        /**
         * Compute the index of the next nonzero kmer.
         *
         * @param i	index of the kmer after which to start searching
         */
        private void findAfter(int i) {
            this.nextKmer = new DnaKmer(i + 1);
            while (this.nextKmer.idx() < size && getBest(this.nextKmer) == Frame.XX)
                nextKmer.increment();
        }

        @Override
        public boolean hasNext() {
            return this.nextKmer.idx() < size;
        }

        @Override
        public DnaKmer next() {
            if (this.nextKmer.idx() >= size)
                throw new NoSuchElementException("Attempt to search past end of Kmer Frame Counter array.");
            // Save the next kmer.
            thisKmer.setIdx(nextKmer.idx());
            // Position after it.
            this.findAfter(this.nextKmer.idx());
            return thisKmer;
        }

    }


    /**
     * Count all of the kmers in a specified genome.
     *
     * @param genome	the genome whose kmers are to be counted
     */
    public void processGenome(Genome genome) {
        // This will hold the kmer inverse.
        DnaKmer invKmer = new DnaKmer();
        // Get the map of location lists.
        Map<String, LocationList> contigMap = LocationList.createGenomeCodingMap(genome);
        // Loop through the contigs from the genome.
        for (Contig contig : genome.getContigs()) {
            // Get the location list for this contig.
            LocationList contigLocs = contigMap.get(contig.getId());
            SequenceDnaKmers kmerProcessor = new SequenceDnaKmers(contig.getSequence());
            while (kmerProcessor.nextKmer()) {
                int pos = kmerProcessor.getPos();
                Frame kmerFrame = contigLocs.computeRegionFrame(pos, pos + DnaKmer.getSize() - 1);
                if (kmerFrame != Frame.XX) {
                    this.increment(kmerProcessor, kmerFrame);
                    invKmer.setRev(kmerProcessor);
                    this.increment(invKmer, kmerFrame.rev());
                }
            }
        }

    }

    /**
     * Erase all the counts so we can start over.
     */
    public void clear() {
        for (int i = 0; i < 7; i++) {
            Arrays.fill(this.countArray[i], (short) 0);
        }
    }

}
