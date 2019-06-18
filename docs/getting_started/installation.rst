Installation
============

User installation
-----------------

FlowCraft is available as a bioconda package, which already comes with
nextflow::

    conda install flowcraft

Alternatively, you can install only FlowCraft, via pip::

    pip install flowcraft

You will also need a container engine (see `Container engine`_ below)

Container engine
----------------

All components of FlowCraft are executed in docker containers, which
means that you'll need to have a container engine installed. The container
engines available are the ones supported by Nextflow:

- `Docker`_,
- `Singularity`_
- Shifter (undocumented)

If you already have any one of these installed, you are good to go. If not,
you'll need to install one. We recommend singularity because it does not
require the processes to run on a separate root daemon.

Singularity
:::::::::::

Singularity is available to download and install `here <http://singularity.lbl.gov/install-linux>`_.
Make sure that you have singularity v2.5.x or higher.
Note that singularity should be installed as root and available on the machine(s) that
will be running the nextflow processes.

.. important::

    Singularity is available as a bioconda package. However, conda installs singularity
    in user space without root privileges, which may prevent singularity images from
    being correctly downloaded. **Therefore it is not recommended that you install
    singularity via bioconda**.

Docker
::::::

Docker can be installed following the instructions on the website:
https://www.docker.com/community-edition#/download.
To run docker as anon-root user, you'll need to following the instructions
on the website: https://docs.docker.com/install/linux/linux-postinstall/#manage-docker-as-a-non-root-user


Developer installation
----------------------

If you are looking to contribute to FlowCraft or simply interested in
tweaking it, clone the github repository and its submodule and then run
setup.py::

    git clone https://github.com/assemblerflow/flowcraft.git
    cd flowcraft
    python3 setup.py install


.. _Docker: https://www.nextflow.io/docs/latest/docker.html
.. _Singularity: https://www.nextflow.io/docs/latest/singularity.html

