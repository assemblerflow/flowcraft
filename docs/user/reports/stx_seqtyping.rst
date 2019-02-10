stx seqtyping
-------------

Table data
^^^^^^^^^^

Typing table:
    - **stx_type_seqtyping**: The stx subtypes found using reads via `seq_typing`.
      Results are presented for `stx1` and `stx2` genes separated by "``:``".

For `stx2` gene, `stx2a`, `stx2c` and `stx2d` variants are grouped together as `stx2acd`
due to the fact that all of these subtypes are the most potent ones to cause HUS and
are difficult to separate from each other by the methods in use right now.

``NT`` (none typeable) is returned when it was not possible to determine a stx subtype.

``NA`` (not available) is returned when ``ecoli_stx_subtyping.py`` was not run.

.. image:: ../resources/reports/stx_seqtyping_table_typing.png
    :align: center
