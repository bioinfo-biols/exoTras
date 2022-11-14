Exosomes recognizing
-------------

Here, we used data in github `tests` directory as an example, and shown how exoTras recognizing exosomes in scRNA-seq datasets.
We have generated test data in `h5ad` format in github, and exoTras also supports `10x_mtx` and `h5` data formats. Noted that the input droplet-gene matrix for exoTras should not be filtered; herein, the matrix should come from the *raw_feature_bc_matrix* directory in Cell Ranger *outs*\.

.. code-block:: python

    import scanpy as sc
    import exoTras
    exoTras.exosomes_recognizer(input_path='./tests', sample_file='./tests/sample_file', out_path='./tests')

The first parameter was the path of directory that contains all samples. Because these test files exists in our `tests` directory, so we used `./tests`. The second parameter was the name of each sample in th directory row by row. If your data format is `10x_mtx`\, exoTras can automatically detect the directory of `sample/outs/raw_feature_bc_matrix/`\.
And **out_path** defines the output of exoTras that is one h5ad file, named *raw_exoTras.h5ad*, with exoTras score and exosome classification in the *obs* for all droplets, and one named *exosomes_exoTras.h5ad* with only exosome-containing droplets.

If samples locate in differet directories, we also supports another way to run exoTras.

.. code-block:: python

    import scanpy as sc
    import exoTras
    exoTras.exosomes_recognizer(sample_file='./tests/sample_file', out_path='./tests')

Here, first parameter was the abosulte path of each sample row by row.

