Pipeline configuration
======================

When a nextflow pipeline is built with assemblerflow, a number of configuration
files are automatically generated in the same directory. They are all imported
at the end of the ``nextflow.config`` file and are sorted by their configuration
role. All configuration files are overwritten if you build another pipeline
in the same directory, with the exception of the ``user.config`` file, which
is meant to be a persistent configuration file.

Parameters
----------

The ``params.config`` file includes all available paramenters for the pipeline
and their respective default values. Most of these parameters already contain
sensible defaults.

Resources
---------

The ``resources.config`` file includes the majority of the directives provided
for each process, including ``cpus`` and ``memory``. You'll note that each
process name has a suffix like ``_1_1``, which is a unique process identifier
composed of ``<lane>_<process_number>``. This ensures that even when the same
component is specified multiple times in a pipeline, you'll still be able to
set directives for each one individually.

Containers
----------

The ``containers.config`` file includes the container directive for each
process in the pipeline. These containers are retrieved from dockerhub, if they
do not exist locally yet. You can change the container string to any other
value, but it should point to an image that exist on dockerhub or locally.

Profiles
--------

The ``profiles.config`` file includes a set of pre-made profiles with all
possible combinations of executors and container engines. You can add new ones
or modify existing one.

User configutations
-------------------

The ``user.config`` file is configuration file that is not overwritten when a
new pipeline is build in the same directory. It can contain any configuration
that is supported by nextflow and will overwrite all other configuration files.