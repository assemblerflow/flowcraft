process_skesa
-------------

Table data
^^^^^^^^^^

Quality control table:
    - **Contigs (skesa)**: Number of assembled contigs.
    - **Assembled BP**: Total number of assembled base pairs.

.. image:: ../resources/reports/assembly_table_skesa.png
    :scale: 80 %
    :align: center

Warnings
^^^^^^^^

Assembly table:
    - When the number of contigs exceeds the threshold of 100 contigs per 1.5Mb.

Fails
^^^^^

Assembly table:
    - When the assembly size if smaller than 80% or larger than 150% of the
      expected genome size.
