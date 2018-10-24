prokka
======


Purpose
-------

This component performs annotations using the annotations available in
`prokka <https://github.com/tseemann/prokka>`_.


Input/Output type
-----------------

- Input type: ``fasta``
- Output type: ``None``

.. note::
    - Although the component doesn't have an output channel it writes the results into the ``publishDir``.


Parameters
----------

- ``centre``: sets the center to which the sequencing center id.
  Default: 'UMMI'.

- ``kingdom``: Selects the annotation mode between Archaea, Bacteria,
  Mitochondria, Viruses. Default: Bacteria).

- ``genus``: Allows user to select a genus name. Default: 'Genus' (same
  as prokka). This also adds the use of the --usegenus flag to prokka.


Published results
-----------------

- ``results/annotation/prokka_<pid>/<sample_id>``: All the outputs from prokka
  will be available in these directories.


Published reports
-----------------

None.


Default directives
------------------

- ``prokka``:
    - ``cpus``: 2
    - ``container``: ummidock/prokka
    - ``version``: 1.12
