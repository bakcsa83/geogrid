package org.giscience.utils.geogrid.cells;

import org.giscience.utils.geogrid.geometry.GeoCoordinates;

/**
 * Grid cell
 *
 * @author Franz-Benjamin Mocnik
 */
public class GridCell implements Comparable<GridCell> {
    private static final double _precision = 1e-9;
    private static final double _precisionPerDefinition = .5e-6;
    private final Integer _resolution;
    private final double _lat;
    private final double _lon;
    private final boolean _isPentagon;
    private Long _id = null;

    public GridCell(int resolution, double lat, double lon, boolean isPentagon) throws Exception {
        if (resolution < 1 || resolution > 23) throw new Exception("resolution must be between 1 and 23");
        this._resolution = resolution;
        if (lat < -90 || lat > 90) throw new Exception("invalid latitude");
        if (lat < -90 + GridCell._precisionPerDefinition || lat > 90 - GridCell._precisionPerDefinition) lon = 0;
        lon %= 360;
        if (lon > 180) lon -= 360;
        else if (lon < -180) lon += 360;
        this._lat = lat;
        this._lon = lon;
        this._isPentagon = isPentagon;
    }

    public GridCell(int resolution, GeoCoordinates c, boolean isPentagon) throws Exception {
        this(resolution, c.getLat(), c.getLon(), isPentagon);
    }

    public Integer getResolution() {
        return this._resolution;
    }

    public Double getLat() {
        return this._lat;
    }

    public Double getLon() {
        return this._lon;
    }

    public Boolean isPentagon() {
        return this._isPentagon;
    }

    /**
     * Returns the ID of the cell. The ID consists of the following elements:
     * <ul>
     *     <li>The leading sign is positive in case of a hexagon, and negative in case of a pentagon.</li>
     *     <li>This first two digits consist of the resolution incremented by 23 in case of negative latitude, by 46 in
     *     case of negative longitude, and by 69 in case of negative latitude and longitude. In case that the
     *     latitude/longitude is strictly less than .5e-6, or that the difference of longitude to 180 or -180 degrees is
     *     strictly less than .5e-6, the respective sign is always regarded as being positive.</li>
     *     <li>The consecutive digits consist of the latitude, with two pre-decimal and six decimal places.</li>
     *     <li>The consecutive digits consist of the longitude, with three pre-decimal and six decimal places. The
     *     longitude is per definition 0 if the latitude differs from -90 or 90 degrees by strictly less than
     *     .5e-6. The longitude is expected to be greater than -180 and strictly less than 180 degrees.</li>
     * </ul>
     *
     * The ID is only valid for resolution smaller less or equal 32.
     *
     * @return ID of the cell
     */
    public Long getID() {
        return this.getID(GridCellIDType.NON_ADAPTIVE);
    }

    public Long getID(GridCellIDType gridCellIDType) {
        if (this._id == null) {
            int numberOfConsecutiveDigits = GridCellMetaData.getInstance().numberOfConsecutiveDigits(this.getResolution(), gridCellIDType);
            double precisionPerDefinition = .5 * Math.pow(10, -numberOfConsecutiveDigits);
            long sgnLat = (this._lat < 0 && Math.abs(this._lat) >= precisionPerDefinition) ? 23 : 0;
            long sgnLon = (this._lon < 0 && Math.abs(this._lon) >= precisionPerDefinition && 180 - Math.abs(this._lon) >= precisionPerDefinition) ? 46 : 0;
            this._id = (this._isPentagon ? -1 : 1) * ((this._resolution.longValue() + sgnLat + sgnLon) * (long) Math.pow(10, 2 * numberOfConsecutiveDigits + 5) + Math.abs(Math.round((this._lat + this._precision) * Math.pow(10, numberOfConsecutiveDigits))) * (long) Math.pow(10, numberOfConsecutiveDigits + 3) + Math.abs(Math.round((this._lon + this._precision) * Math.pow(10, numberOfConsecutiveDigits))));
        }
        return this._id;
    }

    public boolean equals(Object o) {
        return (o instanceof GridCell) && ((GridCell) o).getID().equals(this.getID());
    }

    public int hashCode() {
        return Long.hashCode(this.getID());
    }

    @Override
    public String toString() {
        return String.format("resolution %d lat %f lon %f - %d", this._resolution, this._lat, this._lon, this.getID());
    }

    @Override
    public int compareTo(GridCell o) {
        int d = Integer.compare(this._resolution, o._resolution);
        if (d != 0) return d;
        d = (Math.abs(this._lat - o._lat) < GridCell._precision) ? 0 : Double.compare(this._lat, o._lat);
        if (d != 0) return d;
        return (Math.abs(this._lon - o._lon) < GridCell._precision) ? 0 : Double.compare(this._lon, o._lon);
    }
}
