package org.giscience.utils.geogrid.grids;

import org.giscience.utils.geogrid.cells.GridCell;
import org.giscience.utils.geogrid.cells.GridCellIDType;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

class CellIDAggregatorToFile implements ICellAggregator<CellIDAggregatorToFile> {
    private ArrayList<Long> cells = new ArrayList<>();
    private List<CellIDAggregatorToFile> caList = new ArrayList();
    private final String filename;
    private GridCellIDType gridCellIDType;
    private Integer face = null;
    private int chunk = 0;
    private static final int chunkSize = 10000000;

    public CellIDAggregatorToFile(String filename, GridCellIDType gridCellIDType) {
        this.filename = filename;
        this.gridCellIDType = gridCellIDType;
    }

    @Override
    public ICellAggregator<CellIDAggregatorToFile> cloneEmpty() {
        return new CellIDAggregatorToFile(this.filename, this.gridCellIDType);
    }

    @Override
    public void add(int face, GridCell c) throws IOException {
        this.face = face;
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
        if (this.cells.size() >= this.chunkSize) this.writeChunkToFile();
    }

    @Override
    public void addAll(CellIDAggregatorToFile ca) {
        this.caList.add(ca);
    }

    @Override
    public int size() {
        return 0;
    }

    @Override
    public boolean contains(GridCell c) {
        return false;
    }

    private void writeChunkToFile() throws IOException {
        Collections.sort(this.cells);
        File file = new File(this.filename + ".face" + this.face + "." + this.chunk);
        try (BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(file))) {
            for (Long a : this.cells) {
                bufferedWriter.append(a.toString());
                bufferedWriter.newLine();
            }
        }
        this.chunk++;
        this.cells = new ArrayList<>();
    }

    public void closeFile() throws IOException {
        if (this.cells.size() > 0) this.writeChunkToFile();
        for (CellIDAggregatorToFile ca : this.caList) ca.closeFile();
    }
}
