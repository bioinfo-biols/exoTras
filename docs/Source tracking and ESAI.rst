Source tracking and ESAI
-------------

With the output of `exoTras.exosomes_recognizer`.
Then, we could source these recoginzed exosomes to original cell type and calculate the exosome secretion activity.
```bash
celltype_exosome_number, adata_exo, adata_combined = exoTras.source_tracker(adata_exo, adata_cell, OBSsample='batch', OBScelltype='celltype')

exo_activity_dat = exoTras.ESAI_celltype(adata_exo, adata_cell, OBSsample='batch', OBScelltype='celltype')
```
In the parameters, the two adata variable represent the exosome- and cell- anndata objects. And `batch` and `celltype` are the indexes of sample names and cell types in the anndata of cells. The output of `source_tracker` is in the 'obsm' of 'adata_exo' indexed as 'source'.