Docker containers guidelines
============================

All FlowCraft components require a docker container in order to be executed,
thus if a new component is added, a docker image should be added as well and
uploaded to
.. _docker hub: https://hub.docker.com/ in order to be available to pull in
other machines. Although this can be done in any personal
repository, we recommend that this docker images are added to an already
existing .. _FlowCraft github repository: https://github.com/assemblerflow/docker-imgs
(called here ``Official``) so that docker builds can be automated with github
integration. Also, the centralization of all images will allow other
contributors to easily access and edit these containers instead of forking from
one side to another every time a container needs to be changed/updated.

Official FlowCraft Docker images
--------------------------------

Writing docker images
:::::::::::::::::::::

Official FlowCraft Docker images are available in
.. _this github repository: https://github.com/assemblerflow/docker-imgs .
If you want to add your image to this repository please fork it and make a
Pull Request (PR) with the requested new image or create an issue asking to be
added to the organization as a contributor.


Building docker images
::::::::::::::::::::::

Then, after the image has been added to the FlowCraft
.. _docker-imgs https://github.com/assemblerflow/docker-imgs
github repository, they can be built through
.. _FlowCraft docker hub https://hub.docker.com/u/flowcraft/dashboard/ .

Tag naming
^^^^^^^^^^

Each time a docker image is built using the automated build of docker hub it
should follow this nomenclature: ``version-patch``.
This is used to avoid the override of previous builds for the same images,
allowing for instance users to use different version of the same software using
the same docker image but with different tags.

- ``Version``: Is a string with tree letters like this: ``1.1.1``. Versions should
change every time a new software is added the container.

- ``Patch``: Is a number that follows a ``-`` after the version. Patches should
change every time a change does not affect
the software inside it. For example, updates to database related files required
by some of the software inside the container.

Unofficial FlowCraft Docker images
----------------------------------

Although we **strongly** recommend that all images are stored in FlowCraft
.. _docker-imgs https://github.com/assemblerflow/docker-imgs github repo, it is
not mandatory to do it. Images can be built in another github repo and
also use another docker hub repository to build the images.
However, do make sure that you define it correctly in the directives of the
process as explained in :ref:`DirectivesAnchor`.
