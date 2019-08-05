/**
 *
 */
package org.theseed.genome.kmers;

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
    protected String getLetters() {
        String retVal = getKmerRegion(0);
        return retVal;
    }

    /**
     * @return the letters in the region spanned by this kmer
     *
     * @param offset	if 0, then the first two of each triple are returned; if 1, then the last two
     * 					of each triple are returned
     */
    private String getKmerRegion(int offset) {
        StringBuilder retVal = new StringBuilder(DnaKmer.getSize());
        // Figure out position of the first character pair that won't be in the kmer.
        int n = this.pos + DnaKmer.getSize() / 2 * 3 - 1;
        // Insure we stay in bounds.
        int limit = sequence.length() - 3;
        if (n > limit) n = limit;
        // Loop through getting the letters.
        for (int i = this.pos-1; i < n; i += 3) {
            retVal.append(this.sequence.charAt(i + offset));
            retVal.append(this.sequence.charAt(i + offset + 1));
        }
        return retVal.toString();
    }

    @Override
    public void reverse() {
        // Get the letters at the current position.
        String myLetters = this.getKmerRegion(1);
        this.setIdx(DnaKmer.fromRString(myLetters));
    }

    @Override
    public int regionSize() {
        return DnaKmer.getSize() / 2 * 3;
    }


}
