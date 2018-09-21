process_viral_assembly
----------------------

Table data
^^^^^^^^^^

Quality control table:
    - **Contigs (SPAdes)**: Number of assembled contigs.
    - **Assembled BP (SPAdes)**: Total number of assembled base pairs.
    - **ORFs**: Number of complete ORFs in the assembly.
    - **Contigs (MEGAHIT)**: Number of assembled contigs.
    - **Assembled BP (MEGAHIT)**: Total number of assembled base pairs.


.. image:: ../resources/reports/assembly_table_viral_assembly.png
    :align: center

Fails
^^^^^

Assembly table:
    - When the assembly size if smaller than 80% or larger than 150% of the
      expected genome size.
