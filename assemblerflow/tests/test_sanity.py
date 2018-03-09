import pytest
from contextlib import contextmanager

try:
    import generator.pipeline_parser as ps
    from generator.error_handling import SanityError
except ImportError:
    import assemblerflow.generator.pipeline_parser as ps
    from assemblerflow.generator.error_handling import SanityError


@contextmanager
def not_raises(exception, msg):
    try:
        yield
    except exception:
        raise pytest.fail(msg)


def test_no_brackets():

    pipeline_strs = [
        "A B C | D",
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.brackets_but_no_lanes(p)


def test_no_brackets_pass():

    pipeline_strs = [
        "A B",
        "(A | B)"
    ]

    for p in pipeline_strs:
        with not_raises(SanityError, "pipeline: {}".format(p)):
            ps.brackets_but_no_lanes(p)
