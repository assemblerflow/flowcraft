fastqc
------

Plot data
^^^^^^^^^

- **Base sequence quality**: The average quality score across the read length.

.. image:: ../resources/reports/fastqc_base_sequence_quality.png

- **Sequence quality**: Distribution of the mean sequence quality score.

.. image:: ../resources/reports/fastqc_per_base_sequence_quality.png

- **Base GC content**: Distribution of the GC content of each sequence.

.. image:: ../resources/reports/fastqc_base_gc_content.png

- **Sequence length**: Distribution of the read sequence length.

.. image:: ../resources/reports/fastqc_sequence_length.png

- **Missing data**: Normalized count of missing data across the read length.

.. image:: ../resources/reports/fastqc_missing_data.png


Warnings
^^^^^^^^

The following FastQC categories will issue a warning when they have a ``WARN`` flag:
    - Per base sequence quality.
    - Overrepresented sequences.

The following FastQC categories will issue a warning when do not have a ``PASS`` flag:
    - Per base sequence content.

Fails
^^^^^

The following FastQC categories will issue a fail when they have  a ``FAIL`` flag:
    - Per base sequence quality.
    - Overrepresented sequences.
    - Sequence length distribution.
    - Per sequence GC content.

The following FastQC categories will issue a fail when the do not have a ``PASS`` flag:
    - Per base N content.
    - Adapter content.
