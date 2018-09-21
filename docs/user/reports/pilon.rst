pilon
-----

Table data
^^^^^^^^^^

Quality control table:
    - **Contigs**: Number of assembled contigs.
    - **Assembled BP**: Total number of assembled base pairs.

.. image:: ../resources/reports/assembly_table_skesa.png
    :scale: 80 %
    :align: center

Plot data
^^^^^^^^^

- **Contig size distribution**: Distribution of the size of each assembled contig.

.. image:: ../resources/reports/contig_size_distribution.png

- **Sliding window coverage and GC content**: Provides coverage and GC content
  metrics along the genome using a sliding window approach and two synchronised
  charts.

.. image:: ../resources/reports/sliding_window_amr.png

Warnings
^^^^^^^^

Quality control table:
    - When the enconding and phred score cannot be guessed from the FastQ file(s).

Fails
^^^^^

Quality control table:
    - When the sample has lower estimated coverage than the provided coverage threshold.