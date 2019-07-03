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

}
