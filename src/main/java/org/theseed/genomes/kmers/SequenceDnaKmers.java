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
        this.pos = 1;
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
    public abstract boolean nextKmer();

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

}
