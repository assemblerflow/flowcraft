import pytest

import flowcraft.generator.process_details as pd
import flowcraft.flowcraft as af

from flowcraft.generator.process_collector import collect_process_map
from flowcraft.generator.process_details import COLORS

process_map = collect_process_map()


def test_color_print():

    for c in COLORS:
        pd.colored_print("teste_msg", c)

    assert 1


def test_long_list():

    arguments = af.get_args(["build", "-L"])

    pipeline_string = "fastqc trimmomatic"

    with pytest.raises(SystemExit):
        pd.proc_collector(process_map, arguments, pipeline_string)


def test_short_list():

    arguments = af.get_args(["build", "-l"])

    pipeline_string = "fastqc trimmomatic"

    with pytest.raises(SystemExit):
        pd.proc_collector(process_map, arguments, pipeline_string)
