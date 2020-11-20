package org.giscience.utils.geogrid.grids;

import org.giscience.utils.geogrid.cells.GridCell;

import java.util.HashSet;
import java.util.Set;

class CellAggregator implements ICellAggregator<CellAggregator> {
    private Set<GridCell> cells = new HashSet();

    @Override
    public ICellAggregator<CellAggregator> cloneEmpty() {
        return new CellAggregator();
    }

    @Override
    public void add(int face, GridCell c) {
        this.cells.add(c);
    }

    @Override
    public void addAll(CellAggregator ca) {
        this.cells.addAll(ca.getCells());
    }

    @Override
    public int size() {
        return this.cells.size();
    }

    @Override
    public boolean contains(GridCell c) {
        return this.cells.contains(c);
    }

    public Set<GridCell> getCells() {
        return this.cells;
    }
}
