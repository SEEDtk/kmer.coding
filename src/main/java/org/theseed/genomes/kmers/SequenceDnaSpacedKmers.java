/**
 *
 */
package org.theseed.genomes.kmers;

/**
 * This class handles spaced kmers.  A spaced kmer pulls the first two of three characters out of
 * the contig until it fills the desired size.
 *
 * @author Bruce Parrello
 *
 */
public class SequenceDnaSpacedKmers extends SequenceDnaKmers {

    public SequenceDnaSpacedKmers() {
        super();
    }

    public SequenceDnaSpacedKmers(String sequence) {
        super(sequence);
    }

    @Override
    protected String getLetters(String sequence, int pos) {
        StringBuilder retVal = new StringBuilder(DnaKmer.getSize());
        // Figure out position of the first character pair that won't be in the kmer.
        int n = pos + DnaKmer.getSize() / 2 * 3 - 1;
        int limit = sequence.length() - 2;
        if (n >= limit) n = limit;
        for (int i = pos-1; i < n; i += 3) {
            retVal.append(sequence.charAt(i));
            retVal.append(sequence.charAt(i+1));
        }
        return retVal.toString();
    }


}
