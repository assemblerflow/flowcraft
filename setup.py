import assemblerflow

from setuptools import setup

VERSION = assemblerflow.__version__

with open("README.md") as fh:
    README = fh.read()

setup(
    name="assemblerflow",
    version="{}-1".format(VERSION),
    packages=["assemblerflow",
              "assemblerflow.templates",
              "assemblerflow.generator",
              "assemblerflow.generator.components"],
    package_dir={"assemblerflow": "assemblerflow"},
    package_data={"assemblerflow": ["nextflow.config",
                                    "profiles.config",
                                    "bin/*",
                                    "lib/*",
                                    "generator/templates/*"]},
    data_files=[("", ["LICENSE"])],
    install_requires=[
        "argparse",
        "jinja2",
    ],
    description="A Nextflow pipeline assembler for genomics. Pick your "
                "modules. Assemble them. Run the pipeline.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ODiogoSilva/assemblerflow",
    author="Diogo N Silva",
    author_email="o.diogosilva@gmail.com",
    license="GPL3",
    entry_points={
        "console_scripts": [
            "assemblerflow = assemblerflow.assemblerflow:main"
        ]
    }
)
