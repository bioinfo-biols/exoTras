Functions
---------

The main functions in exoTras are listed below:

.. code-block:: python

    exoTras.exosomes_recognizer(sample_file, out_path, input_path=None, species='Homo', predefine_threads=-2, get_only=False, score_t = None)


.. code-block:: python

    exoTras.source_tracker(adata_exo, adata_cell, OBSsample='batch', OBScelltype='celltype', OBSexo='exo', OBSMpca='X_pca', cellN=10, Xraw = True, normalW=True)


.. code-block:: python

    exoTras.ESAI_celltype(adata_exo, adata_cell, OBSsample='batch', OBScelltype='celltype')

