package org.giscience.utils.geogrid.grids;

import java.util.HashSet;
import java.util.Set;

class ResultCellForBound<T extends ICellAggregator> {
    public final T cellAggregator;
    public final Set<Integer> visitedCells = new HashSet();

    public ResultCellForBound(T cellAggregator) {
        this.cellAggregator = cellAggregator;
    }
}
