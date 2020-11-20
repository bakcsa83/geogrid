package org.giscience.utils.geogrid.grids;

import org.giscience.utils.geogrid.cells.GridCell;

interface ICellAggregator<T> {
    ICellAggregator<T> cloneEmpty();

    void add(int face, GridCell c) throws Exception;

    void addAll(T ca);

    int size();

    boolean contains(GridCell c);
}
