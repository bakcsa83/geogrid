package org.giscience.utils.geogrid.cells;

import org.giscience.utils.geogrid.ISEA3HException;
import org.giscience.utils.geogrid.geometry.GeoCoordinates;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

/**
 * Grid cell
 *
 * @author Franz-Benjamin Mocnik
 */
public class GridCell implements Comparable<GridCell>, Serializable {
    private static final double precision = 1e-9;
    private static final double precisionPerDefinition = 1e-5;
    private final Integer resolution;
    private final double lat;
    private final double lon;
    private final boolean isPentagon;
    private final Map<GridCellIDType, Long> id = new HashMap();

    public GridCell(int resolution, double lat, double lon, boolean isPentagon) throws ISEA3HException {
        if (resolution < 1 || resolution > 22) throw new ISEA3HException("resolution must be between 1 and 22");
        this.resolution = resolution;
        if (lat < -90 || lat > 90) throw new ISEA3HException("invalid latitude");
        if (lat < -90 + GridCell.precisionPerDefinition || lat > 90 - GridCell.precisionPerDefinition) lon = 0;
        lon %= 360;
        if (lon > 180) lon -= 360;
        else if (lon < -180) lon += 360;
        this.lat = lat;
        this.lon = lon;
        this.isPentagon = isPentagon;
    }

    public GridCell(int resolution, GeoCoordinates c, boolean isPentagon){
        this(resolution, c.getLat(), c.getLon(), isPentagon);
    }

    public Integer getResolution() {
        return this.resolution;
    }

    public Double getLat() {
        return this.lat;
    }

    public Double getLon() {
        return this.lon;
    }

    public Boolean isPentagon() {
        return this.isPentagon;
    }

    /**
     * Returns the ID of the cell. The ID consists of the following elements:
     * <ul>
     *     <li>The leading sign is positive in case of a hexagon, and negative in case of a pentagon.</li>
     *     <li>This first two digits consist of the resolution incremented by 22 in case of negative latitude, by 44 in
     *     case of negative longitude, and by 66 in case of negative latitude and longitude. In case that the
     *     latitude/longitude is strictly less than .5e-6, or that the difference of longitude to 180 or -180 degrees is
     *     strictly less than .5e-6, the respective sign is always regarded as being positive.</li>
     *     <li>The consecutive digits consist of the latitude, with two pre-decimal and six decimal places.</li>
     *     <li>The consecutive digits consist of the longitude, with three pre-decimal and six decimal places. The
     *     longitude is per definition 0 if the latitude differs from -90 or 90 degrees by strictly less than
     *     .5e-6. The longitude is expected to be greater than -180 and strictly less than 180 degrees.</li>
     * </ul>
     * <p>
     * If the grid cell type ADAPTIVE_UNIQUE is provided, the values for latitude and longitude are only encoded with
     * the precision required to guarantee uniqueness of the IDs.  In case of ADAPTIVE_1_PERCENT, the values for
     * latitude and longitude are encoded with two more decimal places as in the previous case.
     * <p>
     * The ID is only valid for resolution smaller less or equal 22.
     *
     * @return ID of the cell
     */
    public Long getID() {
        return this.getID(GridCellIDType.NON_ADAPTIVE);
    }

    public Long getID(GridCellIDType gridCellIDType) {
        Long id = this.id.getOrDefault(gridCellIDType, null);
        if (id == null) {
            int numberOfDecimalPlaces = GridCellMetaData.getInstance().numberOfDecimalPlaces(this.getResolution(), gridCellIDType);
            double precisionPerDefinition = .5 * Math.pow(10, -numberOfDecimalPlaces);
            long sgnLat = (this.lat <= -precisionPerDefinition) ? 22 : 0;
            long sgnLon = (this.lon <= -precisionPerDefinition && 180 - Math.abs(this.lon) >= precisionPerDefinition) ? 44 : 0;
            id = (this.isPentagon ? -1 : 1) * ((this.resolution.longValue() + sgnLat + sgnLon) * (long) Math.pow(10, 2 * numberOfDecimalPlaces + 5) + Math.abs(Math.round((this.lat + precision) * Math.pow(10, numberOfDecimalPlaces))) * (long) Math.pow(10, numberOfDecimalPlaces + 3) + Math.abs(Math.round((this.lon + precision) * Math.pow(10, numberOfDecimalPlaces))));
            this.id.put(gridCellIDType, id);
        }
        return id;
    }

    public boolean equals(Object o) {
        return (o instanceof GridCell) && ((GridCell) o).getID().equals(this.getID());
    }

    public int hashCode() {
        return Long.hashCode(this.getID());
    }

    @Override
    public String toString() {
        return String.format("resolution: %d; lat: %f; lon: %f; ID: %d", this.resolution, this.lat, this.lon, this.getID());
    }

    @Override
    public int compareTo(GridCell o) {
        int d = Integer.compare(this.resolution, o.resolution);
        if (d != 0) return d;
        d = (Math.abs(this.lat - o.lat) < GridCell.precision) ? 0 : Double.compare(this.lat, o.lat);
        if (d != 0) return d;
        return (Math.abs(this.lon - o.lon) < GridCell.precision) ? 0 : Double.compare(this.lon, o.lon);
    }
}
