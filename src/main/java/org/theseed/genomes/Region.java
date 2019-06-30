/**
 *
 */
package org.theseed.genomes;

/** This class implements a single region on a contig.  Only the left and right positions (1-based) are given.
 *
 * @author Bruce Parrello
 *
 */
public class Region implements Comparable<Region> {

    // FIELDS
    private int left;
    private int right;

    /**
     * Construct a region from the two positions.
     */
    public Region(int left, int right) {
        this.left = left;
        this.right = right;
    }

    /**
     * @return the left position
     */
    public int getLeft() {
        return left;
    }

    /**
     * @return the right position
     */
    public int getRight() {
        return right;
    }

    /**
     * @param left	the left position to set
     */
    public void setLeft(int left) {
        this.left = left;
    }

    /**
     * @param right the right position to set
     */
    public void setRight(int right) {
        this.right = right;
    }

    /**
     * Shift this region one position to the right.
     */
    public void ShiftRight() {
        this.left++;
        this.right++;
    }

    /**
     * Compare two regions.  The left positions are compared first.  If they are equal, the longer
     * region compares first.
     *
     * @param other	the other region to compare
     */
    @Override
    public int compareTo(Region other) {
        int retVal = this.left - other.left;
        if (retVal == 0) {
            retVal = other.right - this.right;
        }
        return retVal;
    }

}
