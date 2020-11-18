package org.giscience.utils.geogrid.geometry;

/**
 * Cartesian Coordinates of a location on a face of a platonic solid.
 *
 * @author Franz-Benjamin Mocnik
 */
public class FaceCoordinates {
    private final Integer face;
    private final Double x;
    private final Double y;

    public FaceCoordinates(Integer face, Double x, Double y) {
        this.face = face;
        this.x = x;
        this.y = y;
    }

    public Integer getFace() {
        return this.face;
    }

    public Double getX() {
        return this.x;
    }

    public Double getY() {
        return this.y;
    }

    public Double distanceTo(FaceCoordinates c) {
        return Math.sqrt(Math.pow(this.x - c.getX(), 2) + Math.pow(this.y - c.getY(), 2));
    }
    @Override
    public String toString() {
        return String.format("face %d x %f y %f", this.face, this.x, this.y);
    }
}
