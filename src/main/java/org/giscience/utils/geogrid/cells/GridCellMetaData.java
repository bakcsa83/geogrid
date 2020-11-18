package org.giscience.utils.geogrid.cells;

import org.giscience.utils.geogrid.geo.WGS84;
import org.giscience.utils.geogrid.grids.ISEA3H;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Meta data for a grid cell
 *
 * @author Franz-Benjamin Mocnik
 */
public class GridCellMetaData {
    private static final int MAX_NUMBER_OF_DECIMAL_PLACES = 6;
    private static final Map<Integer, Integer> numberOfDecimalPlaces = new ConcurrentHashMap<>();

    private GridCellMetaData() {
    }

    /**
     * Returns the number of decimal places to be used for the ID for a given resolution and a given ID type
     *
     * @param resolution
     * @param gridCellIDType
     * @return number of decimal places
     */
    public static int numberOfDecimalPlaces(int resolution, GridCellIDType gridCellIDType) {
        if (gridCellIDType == GridCellIDType.NON_ADAPTIVE) return MAX_NUMBER_OF_DECIMAL_PLACES;
        int nodp;
        if (numberOfDecimalPlaces.containsKey(resolution)) nodp = numberOfDecimalPlaces.get(resolution);
        else {
            double distBetweenCells = 2 * (new ISEA3H(resolution)).lowerBoundForLengthOfASideOfHexagonalCellOnSphere();
            nodp = (int) Math.ceil(-Math.log10(distBetweenCells / (2 * Math.PI * WGS84.radiusAuthalic / 360)));
            numberOfDecimalPlaces.put(resolution, nodp);
        }
        switch (gridCellIDType) {
            case ADAPTIVE_UNIQUE:
                return Math.max(Math.min(nodp, MAX_NUMBER_OF_DECIMAL_PLACES), 0);
            case ADAPTIVE_1_PERCENT:
                return Math.max(Math.min(nodp + 2, MAX_NUMBER_OF_DECIMAL_PLACES), 0);
            default:
                return MAX_NUMBER_OF_DECIMAL_PLACES;
        }
    }
}
