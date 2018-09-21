import flowcraft

from setuptools import setup

VERSION = flowcraft.__version__

with open("README.md") as fh:
    README = fh.read()

setup(
    name="flowcraft",
    version="{}".format(VERSION),
    packages=["flowcraft",
              "flowcraft.templates",
              "flowcraft.templates.flowcraft_utils",
              "flowcraft.generator",
              "flowcraft.generator.components"],
    package_dir={"flowcraft": "flowcraft"},
    package_data={"flowcraft": ["nextflow.config",
                                "profiles.config",
                                "bin/*",
                                "lib/*",
                                "resources/*",
                                "generator/templates/*"]},
    data_files=[("", ["LICENSE"])],
    install_requires=[
        "pympler",
        "python-dateutil",
        "argparse",
        "jinja2",
        "requests"
    ],
    description="A Nextflow pipeline assembler for genomics. Pick your "
                "modules. Assemble them. Run the pipeline.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/assemblerflow/flowcraft",
    author="Diogo N Silva",
    author_email="o.diogosilva@gmail.com",
    license="GPL3",
    entry_points={
        "console_scripts": [
            "flowcraft = flowcraft.flowcraft:main"
        ]
    }
)
