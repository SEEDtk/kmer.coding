/**
 *
 */
package org.theseed.genomes;

import java.util.TreeSet;

/**
 * @author Bruce Parrello
 *
 * This object maintains a sorted list of non-overlapping locations on a contig.  The
 * locations are inserted by adding features from a genome.  If the feature is on a
 * different contig, it will be rejected.  Features with multiple regions will be
 * marked as invalid.  Overlapping locations will be separated out and marked as invalid.
 *
 * The resulting list can be used to compute coding frames for any kmer location on the
 * contig.
 *
 * Several things are always true about the locations in this list.  They are all on the same
 * contig.  They never overlap.  They are all single-segment.  These assumptions simplify the
 * operations on this object.
 *
 */
public class LocationList {


    // FIELDS
    /** ID of the contig of interest */
    String contigId;
    /** sorted list of locations */
    TreeSet<Location> locations;

    /**
     * Construct a location list for a specified contig.
     *
     * @param contigId	ID of the contig to which this location list is specific
     */
    public LocationList(String contigId) {
        this.contigId = contigId;
        this.locations = new TreeSet<Location>();
    }

    /**
     * Add a new feature to this location list.
     *
     * @param feat	feature to add
     *
     * @return TRUE if the feature was added, FALSE if it belongs to a different contig
     */
    public Boolean addFeature(Feature feat) {
        Boolean retVal = false;
        Location loc = feat.getLocation();
        if (contigId.equals(loc.getContigId())) {
            // Here we are on the same contig, so we can add the location.  Create a single-region
            // copy.
            Location regionLoc = loc.regionOf();
            // Invalidate it if it is segmented.
            if (loc.isSegmented()) {
                regionLoc.invalidate();
            }
            // Now we need to merge it in.  The only tricky part to this is if there is an overlap,
            // we have to create invalid locations for the overlap area. Since there is no overlap
            // inside the list, for another location to overlap with this one, it must be either the
            // floor or ceiling. We remove the adjacent locations, resolve the overlaps, and then
            // add them back.  This list will contain the locations we need to process.
            Location before = this.locations.floor(regionLoc);
            if (before != null && before.getRight() >= regionLoc.getLeft()) {
                locations.remove(before);
                regionLoc = this.ResolveOverlap(before, regionLoc);
            }
            Location after = this.locations.ceiling(regionLoc);
            if (after != null && after.getLeft() <= regionLoc.getRight()) {
                locations.remove(after);
                regionLoc = this.ResolveOverlap(regionLoc, after);
            }
            // Add what's left of the new location.
            locations.add(regionLoc);
            // Denote we've incorporated this location.
            retVal = true;
        }
        return retVal;
    }

    /**
     * Resolve overlaps for a pair of overlapping locations, and add all but the last to the
     * location list.
     *
     * @param loc1	leftmost overlapping location
     * @param loc2	rightmost overlapping location
     *
     * @return a location at the end that has not been incorporated into the list
     */
    private Location ResolveOverlap(Location loc1, Location loc2) {
        Location retVal;
        if (loc1.getRight() > loc2.getRight()) {
            // Here the second location is wholly inside the first.  Mark it as invalid.
            loc2.invalidate();
            // Separate the part of loc1 that precedes loc2.
            Location prefix = loc1.createEmpty();
            prefix.putRegion(loc1.getLeft(), loc2.getLeft() - 1);
            this.locations.add(prefix);
            // Add the new location.
            this.locations.add(loc2);
            // Shorten loc1 and save it as the residual.
            loc1.setLeft(loc2.getRight() + 1);
            retVal = loc1;
        } else if (loc1.getLeft() == loc2.getLeft()) {
            // Here the two locations start at the same place, but the second location
            // extends past the first.  Invalidate the first, and shrink the second.
            loc1.invalidate();
            this.locations.add(loc1);
            loc2.setLeft(loc1.getRight() + 1);
            retVal = loc2;
        } else {
            // Here we have a partial overlap.  We shrink both locations, and create a new
            // one at the end to encompass the suffix.
            Location suffix = loc2.createEmpty();
            suffix.putRegion(loc1.getRight() + 1, loc2.getRight());
            retVal = suffix;
            loc1.setRight(loc2.getLeft() - 1);
            this.locations.add(loc1);
            loc2.setRight(suffix.getLeft() - 1);
            loc2.invalidate();
            this.locations.add(loc2);
        }
        return retVal;
    }

}
