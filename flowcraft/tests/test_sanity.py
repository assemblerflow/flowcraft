import pytest
from contextlib import contextmanager

try:
    import generator.pipeline_parser as ps
    from generator.error_handling import SanityError
except ImportError:
    import flowcraft.generator.pipeline_parser as ps
    from flowcraft.generator.error_handling import SanityError


@contextmanager
def not_raises(exception, msg):
    try:
        yield
    except exception:
        raise pytest.fail(msg)

def test_empty_tasks():
    pipeline_strs = [
        "   ",
        ""
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.empty_tasks(p)


def test_no_brackets_fail():

    pipeline_strs = [
        "A B C | D",
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.brackets_but_no_lanes(p)


def test_number_of_forks_fail():

    pipeline_strs = [
        "A B (( C | D)",
        "A B ( C | D"
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.brackets_insanity_check(p)


def test_lane_char_fail():

    pipeline_strs = [
        "A B (D || E)"
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.lane_char_insanity_check(p)


def test_final_char_fail():

    pipeline_strs = [
        "|",
        "A B |"
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.final_char_insanity_check(p)


def test_fork_no_proc_fail():

    pipeline_strs = [
        "A B (|E)",
        "A B (E|)"
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.fork_procs_insanity_check(p)


def test_double_fork_fail():

    pipeline_strs = [
        "A B (( C | D ) E )"
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.start_proc_insanity_check(p)


def test_close_token_ending_fail():

    pipeline_strs = [
        "A B ( C | D ) E"
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.late_proc_insanity_check(p)


def test_inner_forks_fail():

    pipeline_strs = [
        "A B ( A D )",
    ]

    for p in pipeline_strs:
        with pytest.raises(SanityError):
            ps.inner_fork_insanity_checks(p)


def test_string_pass_all():

    # all these functions listed here don't accept strings with spaces
    pipeline_strs = [
        "A B",
        "(A|B)",
        "A B (C|D)",
        "A B (D|E(F|G))",
        "A B (C|B)",
        "F T(S(P(P|M)|M(P|M(P| M)))|Sp)"
    ]

    for p in pipeline_strs:
        with not_raises(SanityError, "pipeline: {}".format(p)):
            ps.brackets_insanity_check(p)
            ps.lane_char_insanity_check(p)
            ps.brackets_but_no_lanes(p)
            ps.fork_procs_insanity_check(p)
            ps.start_proc_insanity_check(p)
            ps.late_proc_insanity_check(p)


def test_string_spaces_pass_all():

    # this test accepts strings with spaces
    pipeline_strs = [
        "A B",
        "(A | B)",
        "A B ( C | D)",
        "A B (D | E (F | G))",
        "A B ( C | B)",
        # spaces are important for this check
        "F T (S(P(P| M) |M(P|M(P| M)))|Sp)"
    ]

    for p in pipeline_strs:
        with not_raises(SanityError, "pipeline: {}".format(p)):
            ps.inner_fork_insanity_checks(p)


def test_string_pass_all_wrapper():

    pipeline_strs = [
        "A B",
        "(A | B)",
        "A B ( C | D)",
        "A B (D | E (F | G))",
        "A B ( C | B)"
    ]

    for p in pipeline_strs:
        with not_raises(SanityError, "pipeline: {}".format(p)):
            ps.insanity_checks(p)