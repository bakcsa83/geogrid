package org.giscience.utils.geogrid.projections;

import org.giscience.utils.geogrid.generic.Trigonometric;
import org.giscience.utils.geogrid.geometry.FaceCoordinates;
import org.giscience.utils.geogrid.geometry.GeoCoordinates;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * @author Franz-Benjamin Mocnik
 */
public class ISEAProjectionTest {
    private final double precision = 1e-9;
    private final int iterations = 1000000;

    @Test
    public void projectionApartFromPoles() throws Exception {
        ISEAProjection p = new ISEAProjection();
        p.setOrientationSymmetricEquator();
        for (int i = 0; i < this.iterations; i++) {
            GeoCoordinates c = new GeoCoordinates(Math.random() * 179.99 - 89.995, Math.random() * 360);
            FaceCoordinates c2 = p.sphereToIcosahedron(c);
            GeoCoordinates c3 = p.icosahedronToSphere(c2);
            assertTrue(Math.abs(c3.getLat() - c.getLat()) < this.precision);
            assertTrue((Math.abs(c3.getLon() - c.getLon()) % 360) * Trigonometric.cos(c.getLat()) < this.precision);
        }
    }
}
