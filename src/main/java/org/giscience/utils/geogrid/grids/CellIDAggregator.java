package org.giscience.utils.geogrid.grids;

import org.giscience.utils.geogrid.cells.GridCell;
import org.giscience.utils.geogrid.cells.GridCellIDType;

import java.util.HashSet;
import java.util.Set;

class CellIDAggregator implements ICellAggregator<CellIDAggregator> {
    private Set<Long> cells = new HashSet();
    private GridCellIDType gridCellIDType;

    public CellIDAggregator(GridCellIDType gridCellIDType) {
        this.gridCellIDType = gridCellIDType;
    }

    @Override
    public ICellAggregator<CellIDAggregator> cloneEmpty() {
        return new CellIDAggregator(this.gridCellIDType);
    }

    @Override
    public void add(int face, GridCell c) {
        this.cells.add(c.getID(this.gridCellIDType));
    }

    @Override
    public void addAll(CellIDAggregator ca) {
        this.cells.addAll(ca.getCellIDs());
    }

    @Override
    public int size() {
        return this.cells.size();
    }

    @Override
    public boolean contains(GridCell c) {
        return this.cells.contains(c.getID(this.gridCellIDType));
    }

    public Set<Long> getCellIDs() {
        return this.cells;
    }
}
