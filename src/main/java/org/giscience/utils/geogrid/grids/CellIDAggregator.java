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
        switch (this.gridCellIDType) {
            case ADAPTIVE_1_PERCENT:
                this.cells.add(c.getIdAdaptive1Percent());
                break;
            case NON_ADAPTIVE:
                this.cells.add(c.getID());
                break;
            case ADAPTIVE_UNIQUE:
                this.cells.add(c.getIdAdaptiveUnique());
                break;
        }
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
        switch (this.gridCellIDType) {
            case ADAPTIVE_1_PERCENT:
                return this.cells.contains(c.getIdAdaptive1Percent());
            case NON_ADAPTIVE:
                return this.cells.contains(c.getID());
            case ADAPTIVE_UNIQUE:
                return this.cells.contains(c.getIdAdaptiveUnique());
        }
        return false;
    }

    public Set<Long> getCellIDs() {
        return this.cells;
    }
}
