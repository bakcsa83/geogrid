package org.giscience.utils.geogrid.grids;

import org.giscience.utils.geogrid.cells.GridCell;
import org.giscience.utils.geogrid.generic.Tuple;
import org.giscience.utils.geogrid.geo.WGS84;
import org.giscience.utils.geogrid.geometry.FaceCoordinate;
import org.giscience.utils.geogrid.geometry.GeoCoordinates;
import org.giscience.utils.geogrid.projections.ISEAProjection;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * ISEA Aperture 3 Hexagon (ISEA3H) Discrete Global Grid System (DGGS)
 * <p>
 * The ISEA3H grid is constructed by using the icosahedron Snyder equal-area (ISEA) projection to map the surface of the
 * globe to the icosahedron. Thereby, the orientation of the icosahedron is chosen such that the north and the south
 * poles are mapped to the edge midpoints of the icosahedron. The equator is thus mapped symmetrically. A grid
 * (aperture 3) is constructed on the icosahedron, and this grid is mapped back by the inverse projection to the globe.
 * <p>
 * The cells of the grid are identified by the resolution and their center points.
 * <p>
 * The ISEA3H has been proposed by:
 * <p>
 * Kevin Sahr, Denis White, and A. Jon Kimerling: Geodesic Discrete Global Grid Systems. Cartography and Geographic
 * Information Science, 30(2), 121–134, 2003. doi:10._1559/152304003100011090
 *
 * @author Franz-Benjamin Mocnik
 */
public class ISEA3H {
    private static final double precision = 1e-9;
    private final ISEAProjection projection = new ISEAProjection();
    private final int resolution; // resolution
    private final long numberOfHexagonalCells;
    private final int numberOfPentagonalCells = 12;
    private static final double inverseSqrt3 = 1 / Math.sqrt(3); // 1 / \sqrt{3}
    private static final double l0 = ISEAProjection.lengthOfTriangleBase(); // length of the triangle base at resolution 1
    private static final double inverseSqrt3l0 = inverseSqrt3 * l0; // 1 / \sqrt{3} * l0
    private final double l; // length of the triangle base at the given resolution
    private final double l2; // l / 2
    private final double l6; // l / 6
    private final double l23; // l * 2 / 3 - diameter of a hexagonal cell on the icosahedron, in kilometres

    private final double inverseSqrt3l; // 1 / \sqrt{3} * l
    private final double inverseSqrt3l2; // 1 / (2 \sqrt{3}) * l
    private static final double triangleA = l0 / 2.; // l0 / 2 // half base
    private static final double triangleB = inverseSqrt3l0; // 1/\sqrt{3} * l0 // distance center point to tip
    private static final double triangleC = inverseSqrt3l0 / 2.; // 1/(2 \sqrt{3}) * l0 // distance base to center point
    private static final double triangleBCA = (triangleB + triangleC) / triangleA; // (triangleB + triangleC) / triangleA
    private int workerThreadCount = 1;


    private void setDefaultWorkerThreadCount() {
        workerThreadCount = Runtime.getRuntime().availableProcessors() / 2;
        if (workerThreadCount == 0) workerThreadCount = 1;
    }

    public ISEA3H(int resolution) {
        this(resolution, true);
    }

    public ISEA3H(int resolution, boolean rotatedProjection) {
        this(resolution, true, null);
    }

    public ISEA3H(int resolution, boolean rotatedProjection, Integer workerThreadCount) {
        if (rotatedProjection) this.projection.setOrientationSymmetricEquator();
        this.resolution = resolution;
        long numberOfHexagonalCells = 1;
        for (int i = 1; i < this.resolution; i++) numberOfHexagonalCells = 3 * numberOfHexagonalCells + 1;
        this.numberOfHexagonalCells = 20 * numberOfHexagonalCells;

        this.l = Math.pow(inverseSqrt3, this.resolution - 1) * l0;
        this.l2 = this.l / 2.;
        this.l6 = this.l / 6.;
        this.l23 = this.l * 2 / 3.;
        this.inverseSqrt3l = Math.pow(inverseSqrt3, this.resolution) * l0;
        this.inverseSqrt3l2 = this.inverseSqrt3l / 2.;
        if (workerThreadCount != null && workerThreadCount > 0) {
            this.workerThreadCount = workerThreadCount;
        } else {
            setDefaultWorkerThreadCount();
        }

    }

    /**
     * @return diameter of a hexagonal cell on the icosahedron, in kilometres
     */
    public double diameterOfHexagonalCellOnIcosahedron() {
        return this.l23;
    }

    /**
     * @return length of a side of a hexagonal cell on the icosahedron, in kilometres
     */
    public double lengthOfASideOfHexagonalCellOnIcosahedron() {
        return this.l / 3;
    }

    /**
     * @return lower bound for the length of a side of a hexagonal cell on the sphere, in kilometres
     */
    public double lowerBoundForLengthOfASideOfHexagonalCellOnSphere() {
        return ISEAProjection.sphericalDistanceFromCenterToVerticesOnSphere() * WGS84.radiusAuthalic * 2 * Math.PI / (360 * Math.sqrt(Math.pow(3, this.resolution - 1) * 5));
    }

    /**
     * Returns the area of a hexagonal cell. The cells should all have the same area by construction, because the ISEA
     * projection is equal-area.
     *
     * @return area of a hexagonal cell, in square kilometres
     */
    public double areaOfAHexagonalCell() {
        return WGS84.areaOfEarth / (this.numberOfHexagonalCells + 5 / 6. * this.numberOfPentagonalCells);
    }

    /**
     * Returns the area of a pentagonal cell. The cells should all have the same area by construction, because the ISEA
     * projection is equal-area.
     *
     * @return area of a pentagonal cell, in in kilometres
     */
    public double areaOfAPentagonalCell() {
        return 5 / 6. * this.areaOfAHexagonalCell();
    }

    /**
     * @return number of hexagonal cells
     */
    public long numberOfHexagonalCells() {
        return this.numberOfHexagonalCells;
    }

    /**
     * @return number of pentagonal cells
     */
    public int numberOfPentagonalCells() {
        return this.numberOfPentagonalCells;
    }

    /**
     * Returns the grid cell for a given location
     *
     * @param lat latitude
     * @param lon longitude
     * @return corresponding grid cell
     * @throws Exception
     */
    public GridCell cellForLocation(double lat, double lon) throws Exception {
        return this.cellForLocation(new GeoCoordinates(lat, lon));
    }

    /**
     * Returns the grid cell for a given location
     *
     * @param c geographic coordinates
     * @return corresponding grid cell
     * @throws Exception
     */
    public GridCell cellForLocation(GeoCoordinates c) throws Exception {
        FaceCoordinate fc = this.cellForLocation(this.projection.sphereToIcosahedron(c));
        return this.newGridCell(this.projection.icosahedronToSphere(fc), fc);
    }

    /**
     * Returns the coordinates of the center of the corresponding grid cell for given coordinates in the face
     *
     * @param c face coordinates
     * @return face coordinates of the center of the corresponding grid cell
     * @throws Exception
     */
    public FaceCoordinate cellForLocation(FaceCoordinate c) {
        double x = (this.coordinatesNotSwapped()) ? c.getX() : c.getY();
        double y = (this.coordinatesNotSwapped()) ? c.getY() : c.getX();
        int nxCenter = (int) Math.round(x / this.l2);
        double xCenter = nxCenter * this.l2;
        if (Math.abs(x - xCenter) <= this.l6) {
            int nyCenter12 = (int) Math.round(((nxCenter % 2 == 0) ? y : y - this.inverseSqrt3l2) / this.inverseSqrt3l);
            double yCenter12 = nyCenter12 * this.inverseSqrt3l + ((nxCenter % 2 == 0) ? 0 : this.inverseSqrt3l2);
            return this.faceCoordinatesSwapByResolution(c.getFace(), xCenter, yCenter12);
        } else {
            int nyCenter1 = (int) Math.round(y / this.inverseSqrt3l);
            double yCenter1 = nyCenter1 * this.inverseSqrt3l;
            int nyCenter2 = (int) Math.round((y - this.inverseSqrt3l2) / this.inverseSqrt3l);
            double yCenter2 = nyCenter2 * this.inverseSqrt3l + this.inverseSqrt3l2;
            FaceCoordinate cCandidate1 = this.faceCoordinatesSwapByResolution(c.getFace(), xCenter, (nxCenter % 2 == 0) ? yCenter1 : yCenter2);
            FaceCoordinate cCandidate2 = this.faceCoordinatesSwapByResolution(c.getFace(), (x > xCenter) ? xCenter + this.l2 : xCenter - this.l2, (nxCenter % 2 != 0) ? yCenter1 : yCenter2);
            return (c.distance(cCandidate1) < c.distance(cCandidate2)) ? cCandidate1 : cCandidate2;
        }
    }

    /**
     * The following equality holds: cellForLocation = getCoordinatesOfCenter . integerForFaceCoordinates
     * <p>
     * However, the method cellForLocation is defined separately to speed up computations.
     *
     * @param face face to compute the coordinates on
     * @param nx   steps into the direction of the vertex of the hexagon
     * @param ny   steps into the direction of the edge of the hexagon
     * @return coordinates on the face
     */
    private FaceCoordinate getCoordinatesOfCenter(int face, int nx, int ny) {
        double x = nx * this.l2;
        double y = (ny + ((nx % 2 == 0) ? 0 : .5)) * this.inverseSqrt3l;
        return this.faceCoordinatesSwapByResolution(face, x, y);
    }

    private FaceCoordinate getCoordinatesOfCenter(int face, int nx, int ny, int orientation) {
        double x = nx * this.l2;
        double y = (ny - orientation * ((nx % 2 == 0) ? 0 : .5)) * this.inverseSqrt3l;
        return this.faceCoordinatesSwapByResolution(face, x, y);
    }

    private Tuple<Integer, Integer> integerForFaceCoordinates(FaceCoordinate c) {
        double x = (this.coordinatesNotSwapped()) ? c.getX() : c.getY();
        double y = (this.coordinatesNotSwapped()) ? c.getY() : c.getX();
        int nxCenter = (int) Math.round(x / this.l2);
        double xCenter = nxCenter * this.l2;
        if (Math.abs(x - xCenter) <= this.l6)
            return new Tuple(nxCenter, (int) Math.round(((nxCenter % 2 == 0) ? y : y - this.inverseSqrt3l2) / this.inverseSqrt3l));
        else {
            int nyCenter1 = (int) Math.round(y / this.inverseSqrt3l);
            double yCenter1 = nyCenter1 * this.inverseSqrt3l;
            int nyCenter2 = (int) Math.round((y - this.inverseSqrt3l2) / this.inverseSqrt3l);
            double yCenter2 = nyCenter2 * this.inverseSqrt3l + this.inverseSqrt3l2;
            FaceCoordinate cCandidate1 = this.faceCoordinatesSwapByResolution(c.getFace(), xCenter, (nxCenter % 2 == 0) ? yCenter1 : yCenter2);
            FaceCoordinate cCandidate2 = this.faceCoordinatesSwapByResolution(c.getFace(), (x > xCenter) ? xCenter + this.l2 : xCenter - this.l2, (nxCenter % 2 != 0) ? yCenter1 : yCenter2);
            return (c.distance(cCandidate1) < c.distance(cCandidate2)) ? new Tuple(nxCenter, (nxCenter % 2 == 0) ? nyCenter1 : nyCenter2) : new Tuple((x > xCenter) ? nxCenter + 1 : nxCenter - 1, (nxCenter % 2 != 0) ? nyCenter1 : nyCenter2);
        }
    }

    /**
     * Returns the grid cell for the centroid of a given geometry
     *
     * @param g geometry
     * @return corresponding grid cell
     * @throws Exception
     */
    public GridCell cellForCentroid(Geometry g) throws Exception {
        Coordinate c = g.getCentroid().getCoordinate();
        return this.cellForLocation(c.y, c.x);
    }

    /**
     * Returns all cells.
     *
     * @return cells
     */
    public Collection<GridCell> cells() throws Exception {
        return this.cells(new CellAggregator()).getCells();
    }

    /**
     * Returns cells that are inside the bounds, or at least very near.
     * <p>
     * Note that the result should, in fact, include all cells whose center points are inside the given bounds, but also
     * cells nearby. Observe that it can, however, not be guaranteed that all such cells are returned if the longitude
     * range is less than 360 degrees or the latitude range is less than 180 degrees.
     *
     * @param lat0
     * @param lat1
     * @param lon0
     * @param lon1
     * @return cells inside the bounds
     */
    public Collection<GridCell> cellsForBound(double lat0, double lat1, double lon0, double lon1) throws Exception {
        if (lat1 - lat0 >= 180 && lon1 - lon0 >= 360) return this.cells(new CellAggregator()).getCells();
        else return this.cellsForBound(new CellAggregator(), lat0, lat1, lon0, lon1).cellAggregator.getCells();
    }


    private ResultCellForBound cellsForBound(CellAggregator ca, double lat0, double lat1, double lon0, double lon1) throws Exception {
        // compute center of bounding box
        double lat = (lat0 + lat1) / 2.;
        double lon = (lon0 + lon1 + ((lon0 <= lon1) ? 0 : 360)) / 2.;
        if (lon > 360) lon -= 360;
        FaceCoordinate fc = this.projection.sphereToIcosahedron(new GeoCoordinates(lat, lon));
        // compute
        return cellsForBound(new ResultCellForBound(ca), fc, lat0, lat1, lon0, lon1);
    }


    private CellAggregator cells(CellAggregator ca) throws Exception {
        ExecutorService executor = Executors.newFixedThreadPool(this.workerThreadCount);
        List<Future<CellAggregator>> futureList = new ArrayList<>();
        for (int face = 0; face < ISEAProjection.numberOfFaces(); face++) {
            final int f = face;
            futureList.add(executor.submit(new Callable<CellAggregator>() {
                public CellAggregator call() throws Exception {
                    return cellsForFace(new CellAggregator(), f);
                }
            }));
        }

        for (Future<CellAggregator> future : futureList) ca.addAll(future.get().getCells());
        executor.shutdown();
        return ca;
    }

    private ResultCellForBound cellsForBound(ResultCellForBound result, FaceCoordinate fcStart, double lat0, double lat1, double lon0, double lon1) throws Exception {
        // if fcStart is already in result, skip the computation

        if (result.visitedCells.contains(fcStart.getFace())) return result;
        result.visitedCells.add(fcStart.getFace());

        // prepare dNs
        List<Coordinate> dNs = new ArrayList<>();

        dNs.add(new Coordinate(1, 1));
        dNs.add(new Coordinate(-1, 1));
        dNs.add(new Coordinate(1, -1));
        dNs.add(new Coordinate(-1, -1));

        // compute starting coordinates
        Tuple<Integer, Integer> fcn = this.integerForFaceCoordinates(fcStart);
        if (!this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1, fcn._2))) {
            if (this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1 - 1, fcn._2)))
                fcn = new Tuple(fcn._1 - 1, fcn._2);
            else if (this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1 + 1, fcn._2)))
                fcn = new Tuple(fcn._1 + 1, fcn._2);
            else if (this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1, fcn._2 - 1)))
                fcn = new Tuple(fcn._1, fcn._2 - 1);
            else if (this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1, fcn._2 + 1)))
                fcn = new Tuple(fcn._1, fcn._2 + 1);
            else if (this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1 - 1, fcn._2 - 1)))
                fcn = new Tuple(fcn._1 - 1, fcn._2 - 1);
            else if (this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1 + 1, fcn._2 - 1)))
                fcn = new Tuple(fcn._1 + 1, fcn._2 - 1);
            else if (this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1 - 1, fcn._2 + 1)))
                fcn = new Tuple(fcn._1 - 1, fcn._2 + 1);
            else if (this.isCoordinatesInFace(this.getCoordinatesOfCenter(fcStart.getFace(), fcn._1 + 1, fcn._2 + 1)))
                fcn = new Tuple(fcn._1 + 1, fcn._2 + 1);
        }


        for (Coordinate dN : dNs) {
            FaceCoordinate fc;
            GeoCoordinates gc;
            Set<Integer> success = new HashSet();
            Set<Integer> successLast;
            boolean hasFoundInside;
            boolean hasFoundOutsideX = false;
            boolean hasFoundOutsideY = false;
            Map<Integer, FaceCoordinate> faceTodo = new HashMap();
            int nx = fcn._1 - (int) dN.getX() + (((int) dN.getX() < 0) ? (int) dN.getX() : 0);
            while (true) {
                nx += dN.getX();
                int ny = fcn._2 - (int) dN.getY();
                successLast = success;
                Integer maxMinValue = (!successLast.isEmpty()) ? ((dN.getY() >= 0) ? Collections.max(successLast) : Collections.min(successLast)) : null;
                success = new HashSet<>();
                hasFoundInside = false;
                while (true) {
                    ny += dN.getY();
                    fc = this.getCoordinatesOfCenter(fcStart.getFace(), nx, ny);
                    gc = this.projection.icosahedronToSphere(fc);
                    if (this.isInside(gc, lat0, lat1, lon0, lon1)) hasFoundInside = true;
                    else if (hasFoundInside || maxMinValue == null || ((dN.getY() >= 0) ? ny > maxMinValue : ny < maxMinValue))
                        break;
                    else continue;
                    if (this.isCoordinatesInFace(fc)) {
                        success.add(ny);
                        result.cellAggregator.add(this.newGridCell(gc, fc));
                    } else {
                        if ((!hasFoundOutsideX && ny != fcn._2) || (!hasFoundOutsideY && ny == fcn._2)) {
                            if (ny != fcn._2) hasFoundOutsideX = true;
                            else hasFoundOutsideY = true;
                            FaceCoordinate fc2 = this.projection.sphereToIcosahedron(gc);
                            if (faceTodo.containsKey(fc2.getFace())) {
                                faceTodo.put(fc2.getFace(), fc2);
                            } else {
                                int sizeBefore = result.cellAggregator.size();
                                result = this.cellsForBound(result, fc2, lat0, lat1, lon0, lon1);
                                if (result.cellAggregator.size() != sizeBefore) faceTodo.put(fc2.getFace(), null);
                            }
                        }
                        if (!success.isEmpty()) break;
                    }
                }
                if (success.isEmpty()) break;
            }
            for (Map.Entry<Integer, FaceCoordinate> e : faceTodo.entrySet()) {
                if (e.getValue() == null || result.cellAggregator.contains(this.cellForLocation(this.projection.icosahedronToSphere(e.getValue())))) {

                } else {
                    result = this.cellsForBound(result, e.getValue(), lat0, lat1, lon0, lon1);
                }
            }
        }

        return result;
    }

    private CellAggregator cellsForFace(CellAggregator cellAggregator, int face) throws Exception {
        int d = ISEAProjection.faceOrientation(face);
        boolean notSwapped = (this.coordinatesNotSwapped());

        // compute
        FaceCoordinate fc;
        GeoCoordinates gc;
        int nyMax = (int) Math.round((notSwapped ? Math.pow(3, (this.resolution - 1) / 2.) : 2 * Math.pow(3, (this.resolution - 2) / 2.)));
        int nyMin = -(int) (notSwapped ? (nyMax - 1) / 2. : nyMax / 2.);
        int nxMin = 0;
        int nxMax = 0;
        int counter = 0;
        for (int ny = nyMax; ny >= nyMin; ny--) {
            if (notSwapped) {
                if (counter == 1) {
                    nxMax++;
                    nxMin = -nxMax;
                }
                if (counter == 3) {
                    nxMax++;
                    nxMin = -nxMax;
                    counter = 0;
                }
                counter++;
            } else {
                nxMin = (d > 0) ? -(int) Math.floor((nyMax - ny) / 2.) : -(int) Math.ceil((nyMax - ny) / 2.);
                nxMax = (d > 0) ? (int) Math.ceil((nyMax - ny) / 2.) : (int) Math.floor((nyMax - ny) / 2.);
            }
            for (int nx = nxMin; nx <= nxMax; nx++) {
                fc = this.getCoordinatesOfCenter(face, notSwapped ? nx : d * ny, notSwapped ? d * ny : nx, d);
                gc = this.projection.icosahedronToSphere(fc);
                cellAggregator.add(this.newGridCell(gc, fc));
            }
        }

        return cellAggregator;
    }

    private boolean isInside(GeoCoordinates c, double lat0, double lat1, double lon0, double lon1) {
        boolean bLat = lat0 <= c.getLat() && c.getLat() <= lat1;
        boolean bLon = (lon0 <= lon1) ? lon0 <= c.getLon() && c.getLon() <= lon1 : lon0 < c.getLon() || c.getLon() < lon1;
        return bLat && bLon;
    }

    private GridCell newGridCell(GeoCoordinates gc, FaceCoordinate fc) throws Exception {
        int d = this.projection.faceOrientation(fc);
        boolean isPentagon = (Math.abs(Math.abs(fc.getX()) - this.triangleA) < ISEA3H.precision && Math.abs(fc.getY() + d * this.triangleC) < ISEA3H.precision) || (Math.abs(fc.getX()) < ISEA3H.precision && Math.abs(fc.getY() - d * this.triangleB) < ISEA3H.precision);
        return new GridCell(this.resolution, gc, isPentagon);
    }

    /**
     * Returns a buffer for a bounding box of geographic coordinates that needs to be considered in order to ensure
     * that all grid cells, whose center point is within the given bounding box, are contained in the bounding box or
     * the buffer
     *
     * @param lat0
     * @param lat1
     * @param lon0
     * @param lon1
     * @return buffer
     * @throws Exception
     */
    public Double bufferEstimator(double lat0, double lat1, double lon0, double lon1) throws Exception {
        lat0 = Math.max(lat0, -90);
        lat1 = Math.min(lat1, 90);
        if (lon1 - lon0 >= 360) {
            lon0 = -180;
            lon1 = 180;
        }
        FaceCoordinate southWest = this.projection.sphereToIcosahedron(new GeoCoordinates(lat0, lon0));
        FaceCoordinate northEast = this.projection.sphereToIcosahedron(new GeoCoordinates(lat1, lon1));
        GeoCoordinates southWest2 = this.projection.icosahedronToSphere(new FaceCoordinate(southWest.getFace(), southWest.getX() - this.diameterOfHexagonalCellOnIcosahedron() / 2, southWest.getY() - this.diameterOfHexagonalCellOnIcosahedron() / 2));
        GeoCoordinates northEast2 = this.projection.icosahedronToSphere(new FaceCoordinate(northEast.getFace(), northEast.getX() - this.diameterOfHexagonalCellOnIcosahedron() / 2, northEast.getY() - this.diameterOfHexagonalCellOnIcosahedron() / 2));
        List<Double> l = new ArrayList<>();
        l.add(Math.abs(southWest2.getLat() - lat0));
        l.add(Math.abs(southWest2.getLon() - lon0));
        l.add(Math.abs(northEast2.getLat() - lat1));
        l.add(Math.abs(northEast2.getLon() - lon1));
        Optional<Double> result = l.stream().max(Double::compareTo);
        return (result.isPresent()) ? result.get() : null;
    }

    private boolean isCoordinatesInFace(FaceCoordinate c) {
        int d = this.projection.faceOrientation(c);

        // test whether coordinate is left of the triangle, right of the triangle, or below the triangle
        if (d * c.getY() > c.getX() * this.triangleBCA + this.triangleB + ISEA3H.precision) return false;
        if (d * c.getY() > -c.getX() * this.triangleBCA + this.triangleB + ISEA3H.precision) return false;
        if (d * c.getY() < -this.triangleC - ISEA3H.precision) return false;

        return true;
    }

    private FaceCoordinate faceCoordinatesSwapByResolution(int face, double x, double y) {
        return new FaceCoordinate(face, this.coordinatesNotSwapped() ? x : y, this.coordinatesNotSwapped() ? y : x);
    }

    private boolean coordinatesNotSwapped() {
        return this.resolution % 2 != 0;
    }
}
