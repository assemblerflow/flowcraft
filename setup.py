import assemblerflow

from setuptools import setup

VERSION = assemblerflow.__version__

setup(
    name="assemblerflow",
    version=VERSION,
    packages=["assemblerflow",
              "assemblerflow.templates",
              "assemblerflow.generator"],
    package_dir={"assemblerflow": "assemblerflow"},
    package_data={"assemblerflow": ["nextflow.config",
                                    "bin/*",
                                    "lib/*",
                                    "generator/templates/*"]},
    install_requires=[
        "argparse",
        "jinja2",
    ],
    description="Nextflow assembler of bacterial genomic pipelines",
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
