/**
 *
 */
package org.theseed.genomes.kmers.coding;

import org.theseed.genomes.Frame;

/**
 * @author Bruce Parrello
 *
 * This class contains the kmer frame counts for a single kmer.  It is basically an array indexed by frame.
 * We store the numbers as unsigned short integers, which do not exist in java.  To simulate them, the
 * numbers are signed internally, and we return the values as ints with some fancy footwork.  This is all
 * done to insure the program can fit on a laptop.  If we have a big machine available, we can rewrite the
 * fancy stuff to a more natural approach.
 */
public class KmerFrameCounts {

    // FIELDS
    short[] counts;

    /**
     * Construct a blank kmer frame counter.  All the counts will be zero.
     */
    public KmerFrameCounts() {
        this.counts = new short[] { 0, 0, 0, 0, 0, 0, 0 };
    }

    /**
     * @return the kmer count for a specified frame.
     *
     * @param frm	the frame whose count is desired
     */
    public int getCount(Frame frm) {
        return iCount(frm.ordinal());
    }

    /**
     * @return the count for the frame with the specified ordinal
     *
     * @param ordinal	of the frame whose count is desired
     */
    private int iCount(int ordinal) {
        int retVal = 0;
        if (ordinal < 7) {
            retVal = (this.counts[ordinal]) & 0xFFFF;
        }
        return retVal;
    }

    /**
     * Increment the kmer count for a specified frame.
     */
    public void increment(Frame frm) {
        if (frm != Frame.XX) {
            this.counts[frm.ordinal()]++;
        }
    }

    /**
     * @return 	the best frame for this kmer; that is, the frame with the highest count,
     * 			or Frame.XX if the kmer has no instances
     */
    public Frame getBest() {
        int best = 7;
        int bestCount = 0;
        for (int i = 0; i < 7; i++) {
            int count = this.iCount(i);
            if (count > bestCount) {
                bestCount = count;
                best = i;
            }
        }
        return Frame.idxFrame(best);
    }

    /**
     * @return	the fraction of time that the specified frame contains this kmer
     */
    public double getFrac(Frame frm) {
        double retVal = 0;
        int frmCount = this.getCount(frm);
        int total = 0;
        for (int i = 0; i < 7; i++) {
            total += this.iCount(i);
        }
        if (total > 0) {
            retVal = ((double) frmCount) / ((double) total);
        }
        return retVal;
    }

}
