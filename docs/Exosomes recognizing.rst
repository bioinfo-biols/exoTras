Exosomes recognizing
-------------

Use test data in github `tests` as an example.
For a study with multiple scRNA-seq samples, we support two input ways.
.. code-block:: python
    ### first    if samples locate in one directory
    exoTras.exosomes_recognizer(input_path='directory_path', sample_file='sample_name', out_path='output_path', species='Homo')
    ### second   if samples locate in differet directories
    exoTras.exosomes_recognizer(sample_file='sample_path', out_path='output_path', species='Homo')

In the first way, **directory_path** is the path for the directory that contains all sample; and **sample_name** is the file that list each sample in th directory row by row. Every sample should contain directory structure like *sample/outs/raw_feature_bc_matrix/*, and exoTras would automatically recognize for each sample. The second way supports one sample file with abosulte path for each sample, if they are not in one directory.

And **out_path** defines the output of exoTras that is one h5ad file, named *raw_exoTras.h5ad*, with exoTras score and exosome classification in the *obs* for all droplets, and one named *exosomes_exoTras.h5ad* with only exosome-containing droplets.