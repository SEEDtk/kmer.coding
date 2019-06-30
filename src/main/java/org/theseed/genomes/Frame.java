package org.theseed.genomes;

/**
 * Enumeration describing the seven coding frames.  Note the special frame XX indicates an invalid result.
 *
 * @author Bruce Parrello
 *
 */
public enum Frame {
    M3("-3"), M2("-2"), M1("-1"), F0("0"), P1("+1"), P2("+2"), P3("+3"), XX("X");

    /** basic array of all frames */
    private static final Frame[] values = values();

    /** array of all good frames */
    public static final Frame[] all = new Frame[] { M3, M2, M1, F0, P1, P2, P3 };

    /** array size for an array of frame data */
    public static final int nFrames = 7;

    /** reverse complement map */
    private static final Frame[] reverse = new Frame[] { P3, P1, P2, F0, M2, M1, M2 };

    /** modular conversion table for plus strand */
    public static final Frame[] plusFrames = new Frame[] { P3, P1, P2 };

    /** modular conversion table for minus strand */
    public static final Frame[] minusFrames = new Frame[] { M3, M1, M2 };

    /** label of frame */
    String label;

    private Frame(String label) {
        this.label = label;
    }

    /**
     * @return the idx
     */
    public int getIdx() {
        return this.ordinal();
    }

    /**
     * @return the label
     */
    public String toString() {
        return label;
    }

    /**
     * @return the opposite frame (will fail if frame is invalid)
     */
    public Frame rev() {
        return Frame.reverse[this.ordinal()];
    }

    /**
     * @return the frame for a given array index
     *
     * @param idx	array index whose frame object is needed
     */
    public static Frame idxFrame(int idx) {
        return values[idx];
    }

}
