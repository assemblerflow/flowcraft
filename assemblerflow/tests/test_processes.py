import pytest

import assemblerflow.generator.Process as proc


def test_process_init():

    p = proc.Process("init")

    assert [p.template] == ["init"]
