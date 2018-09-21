assembly_mapping
----------------

Plot data
^^^^^^^^^

- **Data loss chart**: Gives a trend of the data loss
  (in total number of base pairs) across components that may filter this data.

.. image:: ../resources/reports/sparkline.png

Warnings
^^^^^^^^

Assembly table:
    - When the number of contigs exceeds the threshold of 100 contigs per 1.5Mb.

Fails
^^^^^

Assembly table:
    - When the assembly size if smaller than 80% or larger than 150% of the
      expected genome size.