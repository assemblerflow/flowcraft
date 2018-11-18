Recipe creation guidelines
===========================

Recipes are pre-made pipeline strings that may be associated with specific
parameters and directives and are used to rapidly build a certain type of
pipeline.

Instead of building a pipeline like::

    -t "integrity_coverage fastqc_trimmomatic fastqc spades pilon"

The user simply can specific a recipe with that pipeline::

    -r assembly

Recipe creation
---------------

The creation of new recipes is a very simple and straightforward process.
You need to create a new file in the ``flowcraft/generator/recipes`` folder
with any name and create a basic class with three attributes::

    try:
        from generator.recipe import Recipe
    except ImportError:
        from flowcraft.generator.recipe import Recipe


    class Innuca(Recipe):

        def __init__(self):
            super().__init__()

            # Recipe name
            self.name = "innuca"

            # Recipe pipeline
            self.pipeline_str = <pipeline string>

            # Recipe parameters and directives
            self.directives = { <directives> }

And that's it! Now there is a new recipe available with the ``innuca`` name and
we can build this pipeline using the option ``-r innuca``.

Name
^^^^

This is the name of the recipe, which is used to make a match with the recipe
name provided by the user via the ``-r`` option.

Pipeline_str
^^^^^^^^^^^^

The pipeline string as if provided via the ``-t`` option.

Directives
^^^^^^^^^^

A dictionary containing the parameters and directives for each process in the
pipeline string. **Setting this attribute is optional and components
that are not specified here will assume their default values**. In general, each
element in this dictionary should have the following format::

    self.directives = {
        "component_name": {
            "params": {
                "paramA": "value"
            },
            "directives": {
                "directiveA": "value"
            }
        }
    }

This will set the provided parameters and directives to the component, but it is
possible to provide only one.

A more concrete example of a real component and directives follows::

    self.pipeline_str = "integrity_coverage fastqc"

    # Set parameters and directives only for integrity_coverage
    # and leave fastqc with the defaults
    self.directives = {
        "integrity_coverage": {
            "params": {
                "minCoverage": 20
            },
            "directives": {
                "memory": "1GB"
            }
        }
    }

Duplicate components
~~~~~~~~~~~~~~~~~~~~

In some cases, the same component may be present multiple times in the pipeline
string of a recipe. In these cases, directives can be assigned to each individual
component by adding a ``#<id>`` suffix to the component::

    self.pipeline_str = "integrity_coverage ( trimmomatic spades#1 | spades#2)"

    self.directives = {
        "spades#1": {
            "directives": {
                "memory": "10GB"
            }
        },
        "spades#2": {
            "directives": {
                "version": "3.7.0"
            }
        }
    }
