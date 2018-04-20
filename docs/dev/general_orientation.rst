General orientation
===================

Codebase structure
------------------

The most important elements of assemblerflow's directory structure are:

- ``generator``:
    - ``components``: Contains the ``Process`` classes for each component
    - ``templates``: Contains the nextflow jinja template files for each component
    - ``engine.py``: The engine of assemblerflow that builds the pipeline
    - ``process.py``: Contains the abstract ``Process`` class that is inherited
    - by all component classes
    - ``pipeline_parser.py``: Functions that parse and check the pipeline string
    - ``recipe.py``: Class responsible for creating recipes
- ``templates``: A git submodule of the `templates`_ repository that contain
  the template scripts for the components.

.. _templates: https://github.com/ODiogoSilva/templates


Code style
----------

- **Style**:  the code base of assemblerflow should adhere (the best it can) to
  the `PEP8`_ style guidelines.
- **Docstrings**: code should be generally well documented following the
  `numpy docstring`_ style.
- **Quality**: there is also an integration with the `codacy`_ service to
  evaluate code quality, which is useful for detecting several coding
  issues that may appear.


Testing
-------

Tests are performed using `pytest`_ and the source files are stored in the
``assemblerflow/tests`` directory. Tests must be executed on the root directory
of the repository

Documentation
-------------

Documentation source files are stored in the ``docs`` directory. The general
configuration file is found in ``docs/conf.py`` and the entry
point to the documentation is ``docs/index.html``.


.. _pytest: https://docs.pytest.org/en/latest/
.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _numpy docstring: https://numpydoc.readthedocs.io/en/latest/format.html
.. _codacy: https://app.codacy.com/app/o.diogosilva/assemblerflow/dashboard