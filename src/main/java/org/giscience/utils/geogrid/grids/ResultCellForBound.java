package org.giscience.utils.geogrid.grids;

import java.util.HashSet;
import java.util.Set;

class ResultCellForBound {
    public final CellAggregator cellAggregator;
    public final Set<Integer> visitedCells = new HashSet();

    public ResultCellForBound(CellAggregator cellAggregator) {
        this.cellAggregator = cellAggregator;
    }
}
