/**
 *
 */
package org.theseed.genome.kmers;

/**
 * This class extracts normal kmers from a contig.  Normal kmers correspond character-for-character
 * to the contig contents.
 *
 * @author Bruce Parrello
 *
 */
public class SequenceDnaNormalKmers extends SequenceDnaKmers {

    public SequenceDnaNormalKmers() {
        super();
    }

    public SequenceDnaNormalKmers(String sequence) {
        super(sequence);
    }

    @Override
    protected String getLetters() {
        int i = this.pos - 1;
        int n = i + DnaKmer.getSize();
        if (n > this.sequence.length()) n = sequence.length();
        return this.sequence.substring(i, n);
    }

    @Override
    public void reverse() {
        this.setIdx(this.rIdx());
    }

    @Override
    public int regionSize() {
        return DnaKmer.getSize();
    }

}
