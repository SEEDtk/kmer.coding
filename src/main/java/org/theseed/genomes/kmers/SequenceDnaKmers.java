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
public class SequenceDnaKmers extends DnaKmer {

    // FIELDS

    /** sequence to iterate through */
    private String sequence;
    /** current position in the sequence */
    private int pos;
    /**
    /**
     * Create a kmer object for extracting all the kmers in a specific DNA sequence.
     *
     * @param sequence	sequence to traverse
     */
    public SequenceDnaKmers(String sequence) {
        super();
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
    public boolean nextKmer() {
        this.pos++;
        int retVal = DnaKmer.fromString(this.sequence, this.pos);
        while (retVal == DnaKmer.NULL) {
            this.pos++;
            retVal = DnaKmer.fromString(this.sequence, this.pos);
        }
        super.setIdx(retVal);
        return (retVal != DnaKmer.EOF);
    }

    /**
     * @return a copy of the kmer at the current position
     */
    public DnaKmer getCopy() {
        return new DnaKmer(this.idx());
    }

}
