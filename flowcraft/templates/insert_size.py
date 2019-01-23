#!/usr/bin/env python3

"""
Purpose
-------

This module is intended to calculate the insert size on sam files
resulting from the mapping the read data back to the assembly.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fasta_file`` : Fasta file paths.
    - e.g.: ``'SampleA.fasta'``
- ``fastq_file``: List with fastq file paths.
    - e.g.: ``'SampleA_1.fastq.gz SampleA_2.fastq.gz'`
- ``plot`` : Boolean to draw insert size plot

Generated output
----------------


Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "09012019"
__template__ = "insert_size-nf"

import os
import sys
import plotly.offline as plot_off
import plotly.graph_objs as graph_obj

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    SAM_STATS = '$sam_stats'
    PLOT = '$plot'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("SAM: {}".format(SAM_STATS))
    logger.debug("PLOT: {}".format(PLOT))


def prepare_insert_size_distribution(samtools_stats):
    """
    Collect the data to produce a distribution plot of the insert sizes

    Parameters
    ----------
    samtools_stats : str
        Path to the samtools stats file

    Returns
    -------
    x : list
        List of x axis values
    y : list
        List of y axis values
    """

    x = []
    y = []
    with open(samtools_stats, 'rt') as reader:
        for line in reader:
            if line.startswith('IS'):
                line = line.split('\t')
                if int(line[2]) > 0:
                    x.append(int(line[1]))
                    y.append(int(line[2]))

    return x, y


def draw_plot(sam_stats, sample_id):

    x_values, y_values = prepare_insert_size_distribution(samtools_stats=sam_stats)

    # Converting absolute counts to frequencies
    y_values = list(map(lambda x: float(x) / sum(y_values), y_values))

    plot_trace = graph_obj.Scatter(name='',
                                   x=x_values,
                                   y=y_values,
                                   mode='lines',
                                   line=dict(color='rgb(0, 0, 0)')
                                   )

    plot_off.plot({"data": [plot_trace],
                   "layout": graph_obj.Layout(title="{} Insert Size Distribution".format(sample_id),
                                              xaxis=dict(title="Insert size"),
                                              yaxis=dict(title="Frequency"))
                   },
                  show_link=False,
                  output_type="file",
                  filename="{}_insert_size_distribution.html".format(sample_id),
                  include_plotlyjs=True,
                  auto_open=False
                  )


@MainWrapper
def main(sample_id, sam_stats, plot):
    """
       Parse Samtools statistics output file and get insert size and standard deviation information

       Parameters
       ----------
       sample_id: str
            Sample name
       sam_stats : str
           Path to the samtools stats file
       plot: Bool
            Boolean to draw or not the insert size plot

       Returns
       -------
       statistics : dict
           Dictionary with the statistics. Keys are statistics description and values are the correspondent values
    """

    statistics = {}

    with open(sam_stats, "rt") as reader:
        counter = 0

        for line in reader:

            if counter <= 2:

                if line.startswith("SN"):
                    line = line.split("\t")

                    if line[1] == "insert size average:":
                        statistics["insert_size"] = round(float(line[2]), 0)
                        counter += 1

                    elif line[1] == "insert size standard deviation:":
                        statistics["insert_size_sd"] = float(line[2])
                        counter += 1
            else:
                break

    with open("{}_insert_size_report.tab".format(sample_id), "wt") as writer:
        writer.write("#" + "\\t".join(sorted(statistics.keys())) + "\\n")
        writer.write("\\t".join([str(statistics[k]) for k in sorted(statistics.keys())]) + "\\n")

    # TODO - .report.json for webApp

    if plot == "True":
        logger.debug("Generating insert size distribution plot.")
        draw_plot(sam_stats, sample_id)



if __name__ == "__main__":

    main(SAMPLE_ID, SAM_STATS, PLOT)

