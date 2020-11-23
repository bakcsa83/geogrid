package org.giscience.utils.geogrid.geometry;

import org.giscience.utils.geogrid.generic.Trigonometric;

/**
 * Geographic coordinates of a location on Earth.
 *
 * @author Franz-Benjamin Mocnik
 */
public class GeoCoordinates implements Comparable<GeoCoordinates> {
    private final Double lat;
    private final Double lon;

    public GeoCoordinates(Double lat, Double lon) throws Exception {
        if (lat < -90 || lat > 90) throw new Exception("invalid latitude");
        lon %= 360;
        if (lon > 180) lon -= 360;
        else if (lon < -180) lon += 360;
        this.lat = lat;
        this.lon = lon;
    }

    public Double getLat() {
        return this.lat;
    }

    public Double getLon() {
        return this.lon;
    }

    public Double distanceTo(GeoCoordinates other) {
        return Trigonometric.acos(Trigonometric.sin(this.getLat()) * Trigonometric.sin(other.getLat()) + Trigonometric.cos(this.getLat()) * Trigonometric.cos(other.getLat()) * Trigonometric.cos(this.getLon() - other.getLon()));
    }

    @Override
    public String toString() {
        return String.format("lat %f lon %f", this.lat, this.lon);
    }

    @Override
    public int compareTo(GeoCoordinates o) {
        int d = Double.compare(this.lat, o.lat);
        if (d != 0) return d;
        return Double.compare(this.lon, o.lon);
    }
}
