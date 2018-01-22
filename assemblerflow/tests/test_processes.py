import pytest

import assemblerflow.generator.Process as proc


def test_process_init():

    p = proc.Process("init", "init")

    assert [p.ptype, p.template] == ["init", "init"]
