import pytest

import assemblerflow.generator.process_details as pd

from assemblerflow.generator.engine import process_map
from assemblerflow.generator.process_details import COLORS


def test_color_print():

    for c in COLORS:
        pd.colored_print("teste_msg", c)

    assert 1


def test_long_list():

    arguments_list = [
        "input_type",
        "output_type",
        "description",
        "dependencies",
        "conflicts"
    ]

    pd.proc_collector(process_map, arguments_list)

    assert 1
