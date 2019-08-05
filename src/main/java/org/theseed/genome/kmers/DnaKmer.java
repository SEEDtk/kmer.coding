/**
 *
 */
package org.theseed.genome.kmers;

/**
 * This class represents a DNA kmer.  All the DNA kmers have a fixed length determined at runtime,
 * and they are represented internally as signed integers.  The special value -1 represents a null
 * kmer.  These particular kmers cannot contain ambiguity characters, and their length cannot be
 * greater than 15.
 *
 * @author Bruce Parrello
 *
 */
public class DnaKmer implements Comparable<DnaKmer> {

    // FIELDS
    /** integer representation of the kmer */
    private int kIdx;

    // RUN-TIME CONSTANTS
    /** number of base pairs in the kmer */
    private static int kmerSize = 15;
    /** null kmer value */
    public static final int NULL = -1;
    /** end-of-string kmer value */
    public static final int EOF = -2;
    /** bit to base converter */
    private static final char[] bitMap = new char[] {'a', 'c', 'g', 't'};


    /**
     * Create a DNA kmer from a DNA string.
     *
     * @param basePairs	DNA sequence to use
     */
    public DnaKmer(String basePairs) {
        this.kIdx = fromString(basePairs, 1);
    }

    /**
     * Create a DNA kmer from a DNA substring.
     *
     * @param basePairs	DNA sequence to use
     * @param pos		position (1-based) in the sequence of the kmer
     */
    public DnaKmer(String basePairs, int pos) {
        this.kIdx = fromString(basePairs, pos);
    }

    /**
     * Create an empty DNA kmer.
     */
    public DnaKmer() {
        this.kIdx = NULL;
    }

    /**
     * Create a DNA kmer from an index number.
     */
    public DnaKmer(int kIdx) {
        this.kIdx = kIdx;
    }

    /**
     * @return the DNA kmer array index.
     */
    public int idx() {
        return this.kIdx;
    }

    /**
     * @return the DNA kmer array index for our reverse compliment.
     */
    public int rIdx() {
        int retVal = 0;
        if (this.kIdx < 0) {
            retVal = NULL;
        } else {
            int reverse = ~this.kIdx;
            for (int i = 0; i < kmerSize; i++) {
                retVal = (retVal << 2) | (reverse & 3);
                reverse >>= 2;
            }
        }
        return retVal;
    }

    /**
     * @return the DNA kmer sequence.
     */
    @Override
    public String toString() {
        return fromKmer(this.kIdx);
    }


    /**
     * Compute the hash code from the kmer's index.
     */
    @Override
    public int hashCode() {
        return this.kIdx;
    }

    /**
     * @return the reverse complement of the DNA kmer sequence.
     */
    public String toRString() {
        int rIdx = this.rIdx();
        return fromKmer(rIdx);
    }

    /**
     * @return the global kmer size.
     */
    public static int getSize() {
        return kmerSize;
    }

    /**
     * Set the global kmer size.
     *
     * @param newSize	new kmer size to use
     */
    public static void setSize(int newSize) {
        kmerSize = newSize;
        if (newSize <= 0) {
            throw new IllegalArgumentException("Invalid kmer size " + newSize + ": cannot be negative.");
        } else if (newSize > 15) {
            throw new IllegalArgumentException("Invalid kmer size " + newSize + ": cannot be greater than 15.");
        }
    }

    /**
     * Store a new value for the kmer index.
     * @param newIdx	new index to store
     */
    public void setIdx(int newIdx) {
        this.kIdx = newIdx;
    }

    // PRIVATE UTILITIES

    /**
     * Compute the string representation of a DNA kmer.  Each base pair is computed from two bits.
     *
     * @param kIdx	the integer index of the kmer
     *
     * @return 	an empty string if the index is negative, else the original kmer from which the index
     * 			was computed
     */
    static private String fromKmer(int kIdx) {
        StringBuilder retVal = new StringBuilder("xxxxxxxxxxxxxxx");
        if (kIdx < 0) {
            retVal.setLength(0);
        } else {
            retVal.setLength(kmerSize);
            for (int i = kmerSize - 1; i >= 0; i--) {
                retVal.setCharAt(i, bitMap[kIdx & 3]);
                kIdx >>= 2;
            }
        }
        return retVal.toString();
    }

    /**
     * Compute the numeric representation of a DNA kmer.  Each base pair is two bits, with 'A' being 0,
     * 'C' being 1, 'G' being 2, and 'T' being 3.
     *
     * @param sequence	the string of base pairs making up the kmer. It is acceptable for it to be
     * 					longer than the kmer size
     * @param pos		the position in the string to start pulling the kmer (1-based)
     *
     * @return the numeric index for the kmer, or -1 if
     */
    static public int fromString(String sequence, int pos) {
        int retVal = 0;
        if ((pos - 1) + kmerSize > sequence.length()) {
            retVal = EOF;
        } else {
            int n = pos + kmerSize - 1;
            for (int i = pos - 1; i < n && retVal >= 0; i++) {
                retVal <<= 2;
                switch (sequence.charAt(i)) {
                case 'a' :
                    break;
                case 'c' :
                    retVal |= 1;
                    break;
                case 'g' :
                    retVal |= 2;
                    break;
                case 't' :
                case 'u' :
                    retVal |= 3;
                    break;
                default :
                    retVal = NULL;
                }
            }
        }
        return retVal;
    }

    /**
     * Compute the numeric representation of the reverse compliment of a DNA kmer.
     * This is more limited than fromString, since it does not allow a starting offset
     * and the string must be the proper size.
     *
     * @param sequence	the string of base pairs making up the kmer.
     *
     * @return the numeric index for the kmer, or -1 if
     */
    static public int fromRString(String sequence) {
        int retVal = 0;
        if (kmerSize > sequence.length()) {
            retVal = EOF;
        } else {
            for (int i = kmerSize - 1; i >= 0 && retVal >= 0; i--) {
                retVal <<= 2;
                switch (sequence.charAt(i)) {
                case 'a' :
                    retVal |= 3;
                    break;
                case 'c' :
                    retVal |= 2;
                    break;
                case 'g' :
                    retVal |= 1;
                    break;
                case 't' :
                case 'u' :
                    break;
                default :
                    retVal = NULL;
                }
            }
        }
        return retVal;
    }

    /**
     * @return the number of possible kmers
     */
    public static int maxKmers() {
        return 4 << (2 * (kmerSize - 1));
    }

    /**
     * Increment this kmer to get the next sequential kmer by index.
     */
    public void increment() {
        this.kIdx++;
    }

    @Override
    public int compareTo(DnaKmer arg0) {
        // Insure EOF compares last.
        int retVal;
        if (this.kIdx == DnaKmer.EOF) {
            if (arg0.kIdx == DnaKmer.EOF) {
                retVal = 0;
            } else {
                retVal = 1;
            }
        } else if (arg0.kIdx == DnaKmer.EOF) {
            retVal = -1;
        } else {
            retVal = (this.kIdx - arg0.kIdx);
        }
        return retVal;
    }

    @Override
    public boolean equals(Object arg0) {
        return (this.kIdx == ((DnaKmer) arg0).kIdx);
    }

    /**
     * @return TRUE if this kmer is the other kmer's reverse complement.
     */
    public boolean isRev(DnaKmer other) {
        return (this.kIdx == other.rIdx());
    }


}
