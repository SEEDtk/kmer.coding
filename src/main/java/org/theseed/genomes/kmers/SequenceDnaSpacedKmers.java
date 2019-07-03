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
    public boolean nextKmer() {
        // TODO get next spaced kmer
        return false;
    }

}
