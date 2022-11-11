### Basic Usage
Input for exoTras is a cell-by-gene matrix. In the case of scRNA-seq dataset using 10X Genomics, we used `raw_feature_bc_matrix` directory generated by Cell Ranger as input. Output of exoTras consists of the score of exosome signals and classification for each droplet. Such exosome information will be used for downstream analysis and as basis for the construction of the exosome secretion activity index (ESAI) and source tracking functions.

We implemented four functions for exosomes recognizing and functional analyses. In exoTras, `exosomes_recognizer` recognizes exosome-containing droplets in the raw scRNA-seq data; `source_tracker` and `ESAI_celltype` trace these droplets to their original cell type and estimited corresponding exosome secretion activity; and `cellfree_simulator` simulates transcriptional profile of cell free droplets in scRNA-seq.