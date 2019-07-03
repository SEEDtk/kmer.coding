/**
 *
 */
package org.theseed.genomes.kmers;

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
    protected String getLetters(String sequence, int pos) {
        int i = pos - 1;
        int n = i + DnaKmer.getSize();
        if (n > sequence.length()) n = sequence.length();
        return sequence.substring(i, n);
    }

}
