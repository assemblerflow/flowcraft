check_coverage
--------------

Table data
^^^^^^^^^^

Quality control table:
    - **Coverage**: Estimated coverage based on the number of base pairs and the expected
      genome size.

.. image:: ../resources/reports/quality_control_table.png
    :align: center

Warnings
^^^^^^^^

Quality control table:
    - When the enconding and phred score cannot be guessed from the FastQ file(s).

Fails
^^^^^

Quality control table:
    - When the sample has lower estimated coverage than the provided coverage threshold.