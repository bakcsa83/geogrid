package org.giscience.utils.geogrid.grids;

import org.giscience.utils.geogrid.cells.GridCell;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public class CellAggregator {
    private Set<GridCell> cells = new HashSet();

    public CellAggregator() {
    }

    void add(GridCell cell) {
        cells.add(cell);
    }

    void addAll(Collection<GridCell> gridCells) {
        cells.addAll(gridCells);
    }

    public int size() {
        return cells.size();
    }

    public boolean contains(GridCell gridCell) {
        return cells.contains(gridCell);
    }

    public Collection<GridCell> getCells() {
        return cells;
    }

}
