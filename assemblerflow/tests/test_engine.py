import pytest

import assemblerflow.generator.engine as eg

from assemblerflow.generator.engine import process_map


@pytest.fixture
def single_con():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc", "lane": 1}}
           ]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def single_fork():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "spades", "lane": 2}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "skesa", "lane": 3}},
           ]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def raw_forks():

    con = [{"input": {"process": "__init__", "lane": 0},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc", "lane": 1}},
           {"input": {"process": "__init__", "lane": 0},
            "output": {"process": "patho_typing", "lane": 2}},
           {"input": {"process": "__init__", "lane": 0},
            "output": {"process": "seq_typing", "lane": 3}}]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def multi_forks():

    con = [{"input": {"process": "__init__", "lane": 0},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "__init__", "lane": 0},
            "output": {"process": "seq_typing", "lane": 2}},
           {"input": {"process": "__init__", "lane": 0},
            "output": {"process": "trimmomatic", "lane": 3}},
           {"input": {"process": "trimmomatic", "lane": 3},
            "output": {"process": "check_coverage", "lane": 3}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "spades", "lane": 4}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "skesa", "lane": 5}},
           {"input": {"process": "check_coverage", "lane": 3},
            "output": {"process": "spades", "lane": 6}},
           {"input": {"process": "check_coverage", "lane": 3},
            "output": {"process": "skesa", "lane": 7}}]

    return eg.NextflowGenerator(con, "teste.nf")


def test_simple_init():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "{}", "lane": 1}}]

    for p in process_map:

        con[0]["output"]["process"] = p
        nf = eg.NextflowGenerator(con, "teste.nf")

        assert [len(nf.processes), nf.processes[1].template] == \
            [2, p]


def test_connections_single_process_channels(single_con):

    template = "integrity_coverage"

    p = single_con.processes[1]

    assert [p.input_channel, p.output_channel] == \
        ["{}_in_1_0".format(template), "{}_out_1_0".format(template)]


def test_connections_invalid():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "spades", "lane": 1}},
           {"input": {"process": "spades", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}}
           ]

    with pytest.raises(SystemExit):
        eg.NextflowGenerator(con, "teste.nf")


def test_connections_nofork(single_con):

    assert single_con._fork_tree == {}


def test_connections_singlefork(single_fork):

    assert single_fork._fork_tree == {1: [2, 3]}


def test_connections_rawfork(raw_forks):

    assert raw_forks._fork_tree == {0: [1, 2, 3]}


def test_connections_multiforks(multi_forks):

    assert multi_forks._fork_tree == {0: [1, 2, 3], 1: [4, 5], 3: [6, 7]}


def test_connections_no_fork_channel_update(single_con):

    p = single_con.processes[1]

    assert p.forks == []


def test_connections_fork_channel_update(single_fork):

    p = single_fork.processes[1]

    assert p.forks != []


def test_connections_channel_update(single_con):

    p1 = single_con.processes[1]
    p2 = single_con.processes[2]

    assert p1.output_channel == p2.input_channel


def test_connections_channel_update_wfork(single_fork):

    p1 = single_fork.processes[1]
    p2 = single_fork.processes[2]
    p3 = single_fork.processes[3]

    assert [p1.main_forks[1], p1.main_forks[2]] == \
           [p2.input_channel, p3.input_channel]
