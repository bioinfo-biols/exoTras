Source tracking and ESAI
-------------

With the output of ``exoTras.exosomes_recognizer`` in :doc:`/Exosomes recognizing` and cell matrix with cell type, exoTras can track each exosome to original cell type and calculate exosome secretion activity index (ESAI).

.. code-block:: python

    celltype_exosome_number, adata_exo, adata_combined = exoTras.source_tracker(adata_exo, adata_cell, Xraw = False)
    exo_activity_dat = exoTras.ESAI_celltype(adata_exo, adata_cell)


The first two parameters represent the exosome- and cell- anndata objects. The last parameter means that do not use the raw object in the ``adata_cell``\. Here, exoTras automatically used information of ``celltype`` in the ``obs`` of ``adata_cell``. We can change the index with the parameters ``OBScelltype``\.

The original cell type for each was listed in the ``obsm`` of ``adata_exo`` indexed as ``source``.
And ``exo_activity_dat`` includes values of ESAI for each cell type as follows:

+------------+------------+-------------+------------+-----------+
|      0     |  celltype0 |  celltype1  |  celltype2 | celltype3 |
+============+============+=============+============+===========+
| test1_h5ad |    4.25    | 0.090909091 |      0     |     0     |
+------------+------------+-------------+------------+-----------+
| test2_h5ad |      0     |    0.0625   |     1.5    |   1.9375  |
+------------+------------+-------------+------------+-----------+
