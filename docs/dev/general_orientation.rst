General orientation
===================

Codebase structure
------------------

The most important elements of assemblerflow's directory structure are:

- ``generator``:
    - ``components``: Contains the ``Process`` classes for each component
    - ``templates``: Contains the nextflow template files for each component
    - ``engine.py``: The engine of assemblerflow that builds the pipeline
    - ``process.py``: Contains the abstract ``Process`` class that is inherited
    - by all component classes
    - ``pipeline_parser.py``: Functions that parse and check the pipeline string
    - ``recipe.py``: Class responsible for creating recipes
- ``templates``: A git submodule of the `templates`_ repository that contain
  the template scripts for the components.

.. _templates: https://github.com/ODiogoSilva/templates