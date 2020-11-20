package org.giscience.utils.geogrid.geometry;

import org.locationtech.jts.geom.Coordinate;

/**
 * Cartesian Coordinates of a location on a face of a platonic solid.
 *
 * @author Franz-Benjamin Mocnik
 */
public class FaceCoordinate extends Coordinate {
    private final Integer face;

    public FaceCoordinate(Integer face, Double x, Double y) {
        this.face = face;
        this.x = x;
        this.y = y;
    }

    public Integer getFace() {
        return this.face;
    }

    @Override
    public String toString() {
        return String.format("face %d x %f y %f", this.face, this.x, this.y);
    }
}
