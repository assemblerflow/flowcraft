import pytest
import os

import flowcraft.generator.utils as utils
from flowcraft.generator.error_handling import LogError


def test_empty_log():
    with pytest.raises(LogError):
        utils.get_nextflow_filepath(
            os.path.join(os.getcwd(), "flowcraft/tests/broadcast_tests/empty_log.txt"))


def test_no_path_in_log():
    with pytest.raises(LogError):
        utils.get_nextflow_filepath(
            os.path.join(os.getcwd(), "flowcraft/tests/broadcast_tests/log_without_command.txt"))


def test_path_in_log():
    filepath = utils.get_nextflow_filepath(
        os.path.join(os.getcwd(), "flowcraft/tests/broadcast_tests/log_with_command.txt"))

    assert filepath != ""


def test_regex_in_log():
    filepath = utils.get_nextflow_filepath(
        os.path.join(os.getcwd(), "flowcraft/tests/broadcast_tests/log_with_command_regex.txt"))

    assert filepath != ""
