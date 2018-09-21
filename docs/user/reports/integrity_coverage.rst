integrity_coverage
------------------

Table data
^^^^^^^^^^

Quality control table:
    - **Raw BP**: Number of raw base pairs from the FastQ file(s).
    - **Reads**: Number of reads in the FastQ file(s)
    - **Coverage**: Estimated coverage based on the number of base pairs and the expected
      genome size.

.. image:: ../resources/reports/quality_control_table.png
    :align: center

Plot data
^^^^^^^^^

- **Data loss chart**: Gives a trend of the data loss
  (in total number of base pairs) across components that may filter this data.

.. image:: ../resources/reports/sparkline.png

Warnings
^^^^^^^^

Quality control table:
    - When the enconding and phred score cannot be guessed from the FastQ file(s).

Fails
^^^^^

Quality control table:
    - When the sample has lower estimated coverage than the provided coverage threshold.