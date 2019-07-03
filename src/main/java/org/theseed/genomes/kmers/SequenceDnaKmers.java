/**
 *
 */
package org.theseed.genomes.kmers;

/**
 * Class for looping through the kmers in a DNA sequence.
 *
 * @author Bruce Parrello
 *
 */
public abstract class SequenceDnaKmers extends DnaKmer {

    // FIELDS

    /** sequence to iterate through */
    protected String sequence;
    /** current position in the sequence */
    protected int pos;

    // Create a blank sequence kmer traversal object.
    protected SequenceDnaKmers() { }


    /**
     * Create a kmer object for extracting all the kmers in a specific DNA sequence.
     *
     * @param sequence	sequence to traverse
     */
    public SequenceDnaKmers(String sequence) {
        super();
        init(sequence);
    }


    /**
     * Initialize this object to traverse the specified sequence.
     *
     * @param sequence	sequence to traverse
     */
    protected void init(String sequence) {
        // Save the sequence.
        this.sequence = sequence;
        // Denote we haven't started.
        this.pos = 0;
    }


    /**
     * @return the current position in the sequence
     */
    public int getPos() {
        return this.pos;
    }

    /**
     * Get the next valid kmer in the sequence.
     *
     * @return TRUE if successful, FALSE if we are at end-of-sequence
     */
    public boolean nextKmer() {
        int retVal = DnaKmer.NULL;
        while (retVal == DnaKmer.NULL) {
            this.pos++;
            String letters = this.getLetters();
            retVal = DnaKmer.fromString(letters, 1);
        }
        super.setIdx(retVal);
        return (retVal != DnaKmer.EOF);
    }

    /**
     * @return the sequence of letters at the current position used to build a kmer of this type.
     */
    protected abstract String getLetters();

    /**
     * @return a copy of the kmer at the current position
     */
    public DnaKmer getCopy() {
        return new DnaKmer(this.idx());
    }

    /**
     * Create a sequence kmer-traversal object for the specified kmer type.
     *
     * @param sequenceClass	class representing the type of kmer to use
     * @param sequence		sequence to traverse
     *
     * @return	an object of the proper subclass type to traverse the kmers in the sequence
     */
    public static SequenceDnaKmers build(Class<? extends SequenceDnaKmers> sequenceClass, String sequence) {
        SequenceDnaKmers retVal;
        try {
            retVal = sequenceClass.newInstance();
            retVal.init(sequence);
        } catch (Exception e) {
            throw new RuntimeException("Error creating kmer-traversal object.", e);
        }
        return retVal;
    }


    /**
     * Set the embedded DNA kmer to the reverse complement of what's at the current position.
     */
    public abstract void reverse();

    /**
     * @return the substring at the current position in the DNA string
     *
     * @param len	length of substring to return
     *
     * NOTE this is usually used for debugging.
     */
    public String stringAtPos(int len) {
        int start = this.pos - 1;
        int end = start + len;
        if (end > this.sequence.length()) end = this.sequence.length();
        return this.sequence.substring(start, end);
    }

    /**
     * @return the size of the region covered by a kmer for this processor
     */

    public abstract int regionSize();

}
