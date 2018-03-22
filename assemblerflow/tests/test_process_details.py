import pytest

import assemblerflow.generator.process_details as pd
import assemblerflow.assemblerflow as af

from assemblerflow.generator.engine import process_map
from assemblerflow.generator.process_details import COLORS


def test_color_print():

    for c in COLORS:
        pd.colored_print("teste_msg", c)

    assert 1


def test_long_list():

    arguments = af.get_args(["-L"])

    with pytest.raises(SystemExit):
        pd.proc_collector(process_map, arguments)


def test_short_list():

    arguments = af.get_args(["-l"])

    with pytest.raises(SystemExit):
        pd.proc_collector(process_map, arguments)
        