package org.giscience.utils.geogrid.grids;

import org.giscience.utils.geogrid.cells.GridCell;
import org.giscience.utils.geogrid.geo.WGS84;
import org.giscience.utils.geogrid.geometry.FaceCoordinate;
import org.junit.Assert;
import org.junit.Test;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;

import java.util.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author Franz-Benjamin Mocnik
 */
public class ISEA3HTest {
    private final double _precision = 1e-9;
    private final double _precision2 = 1e-15;
    private final int _iterations = 1000000;

    @Test
    public void numberOfCellsByNumber() {
        this.numberOfCellsByNumber(false);
    }

    @Test
    public void numberOfCellsByNumberRotated() {
        this.numberOfCellsByNumber(true);
    }

    public void numberOfCellsByNumber(boolean rotatedProjection) {
        this._numberOfCellsByNumber(1, 20, rotatedProjection);
        this._numberOfCellsByNumber(2, 80, rotatedProjection);
        this._numberOfCellsByNumber(3, 260, rotatedProjection);
        this._numberOfCellsByNumber(4, 800, rotatedProjection);
        this._numberOfCellsByNumber(15, 143489060, rotatedProjection);
        this._numberOfCellsByNumber(16, 430467200, rotatedProjection);
    }

    public void _numberOfCellsByNumber(int resolution, int numberOfHexagonalCells, boolean rotatedProjection) {
        ISEA3H grid = new ISEA3H(resolution, rotatedProjection);
        assertEquals(grid.numberOfHexagonalCells(), numberOfHexagonalCells);
        assertEquals(grid.numberOfPentagonalCells(), 12);
    }

    @Test
    public void numberOfCellsByArea() {
        this.numberOfCellsByArea(false);
    }

    @Test
    public void numberOfCellsByAreaRotated() {
        this.numberOfCellsByArea(true);
    }

    public void numberOfCellsByArea(boolean rotatedProjection) {
        for (int r = 1; r < 19; r++) this.numberOfCellsByArea(r, rotatedProjection);
    }

    public void numberOfCellsByArea(int resolution, boolean rotatedProjection) {
        ISEA3H grid = new ISEA3H(resolution, rotatedProjection);
        assertTrue(grid.numberOfHexagonalCells() * grid.areaOfAHexagonalCell() + grid.numberOfPentagonalCells() * grid.areaOfAPentagonalCell() - WGS84.areaOfEarth < WGS84.areaOfEarth * this._precision2);
    }

    @Test
    public void numberOfCellsByGrid() throws Exception {
        this.numberOfCellsByGrid(false);
    }

    @Test
    public void numberOfCellsByGridRotated() throws Exception {
        this.numberOfCellsByGrid(true);
    }

    public void numberOfCellsByGrid(boolean rotatedProjection) throws Exception {
        for (int r = 1; r <= 12; r++) this.numberOfCellsByGrid(r, rotatedProjection);
    }

    public void numberOfCellsByGrid(int resolution, boolean rotatedProjection) throws Exception {
        ISEA3H grid = new ISEA3H(resolution, rotatedProjection);
        long cells = grid.cells().size();
        long cellsExpected = grid.numberOfHexagonalCells() + grid.numberOfPentagonalCells();
        System.out.format("resolution %d - %d - %d\n", resolution, cells, cellsExpected);
        assertEquals(cells, cellsExpected);
    }

    @Test
    public void pointsInGridCells() throws Exception {
        this._pointsInGridCells(1);
        this._pointsInGridCells(12);
        this._pointsInGridCells(13);
        this._pointsInGridCells(16);
    }

    private void _pointsInGridCells(int resolution) throws Exception {
        ISEA3H grid = new ISEA3H(resolution);
        for (int i = 0; i < this._iterations; i++) {
            FaceCoordinate c = new FaceCoordinate(1, Math.random() * 100 - 50, Math.random() * 100 - 50);
            assertTrue(c.distance(grid.cellForLocation(c)) <= grid.diameterOfHexagonalCellOnIcosahedron() / 2. + this._precision);
        }
    }

    @Test
    public void cellTest() throws Exception {
        long start = System.nanoTime();
        ISEA3H g = new ISEA3H(8);
        Collection<GridCell> cells = g.cells();
        Map<Long, GridCell> cellMap = new HashMap<>();

        Iterator<GridCell> i = cells.iterator();
        while (i.hasNext()) {
            GridCell c = i.next();
            cellMap.put(c.getID(), c);
        }

        GridCell testCell = cellMap.get(3018730356082759079L);
        Assert.assertNotNull(testCell);
        Assert.assertEquals(testCell.getLat(), -18.730356348906373, 0.0000000000000001);
        Assert.assertEquals(testCell.getLon(), 82.7590791756899, 0.0000000000000001);

        testCell = cellMap.get(3006205244125280036L);
        Assert.assertNotNull(testCell);
        Assert.assertEquals(testCell.getLat(), -6.205244105979422, 0.0000000000000001);
        Assert.assertEquals(testCell.getLon(), 125.28003578854273, 0.0000000000000001);
    }

    @Test
    public void cellForCoordinateTest() throws Exception {
        ISEA3H g = new ISEA3H(17);
        GridCell cell = g.cellForLocation(45, 150);
        Assert.assertEquals(1744997589150003607L,cell.getID().longValue());
    }

    @Test
    public void coordinateTest() throws Exception {
        int testPointCount = 2000;
        double lonStep = 350.0 / (double) testPointCount;
        double latStep = 170.0 / (double) testPointCount;
        Coordinate startPoint = new Coordinate(-179, -89);
        ISEA3H g = new ISEA3H(12);

        Point[] testPoints = new Point[testPointCount];
        Set<GridCell> cells = new HashSet<>();
        Map<Long, List<Point>> pointToCellMap = new HashMap<>();

        GeometryFactory gf = new GeometryFactory();
        for (int i = 0; i < testPointCount; i++) {
            testPoints[i] = gf.createPoint(new Coordinate(startPoint.getX() + (lonStep * (i + 1)), startPoint.getY() + (latStep * (i + 1))));
            GridCell cell = g.cellForLocation(testPoints[i].getY(), testPoints[i].getX());
            cells.add(cell);
            List<Point> pointList = pointToCellMap.get(cell.getID());
            if (pointList == null) {
                pointList = new ArrayList<>();
                pointList.add(testPoints[i]);
                pointToCellMap.put(cell.getID(), pointList);
            } else {
                pointList.add(testPoints[i]);
            }
        }

        Iterator<GridCell> cellsIterator = cells.iterator();
        Point[] cellCoordinates = new Point[cells.size()];
        int counter = 0;
        while (cellsIterator.hasNext()) {
            GridCell cell = cellsIterator.next();
            cellCoordinates[counter++] = gf.createPoint(new Coordinate(cell.getLon(), cell.getLat()));
        }

        Assert.assertEquals(1991, cellCoordinates.length);

        int mappedPoints=0;
        for(List<Point> list:pointToCellMap.values()){
            mappedPoints+=list.size();
        }
        Assert.assertEquals(testPointCount,mappedPoints);
    }

}
