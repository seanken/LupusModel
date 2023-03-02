# LupusModel

This directory contains code for the manuscript "".

The CellRanger directory contains the script used to run CellRanger on our samples.

The R directory contains the scripts for downstream analysis. The main scripts are:

**LoadData.R**: Code used to load CellRanger output and downsample.

**Load.set2.from.mat.R**: Code used to cluster, ID cell types, etc, to get final Seurat object.

**RunDownstream.set2.R**: Code used for downstream analysis (DE, enrichment, etc).

In addition, the Tools subdirectory contains many functions used by the main scripts.
