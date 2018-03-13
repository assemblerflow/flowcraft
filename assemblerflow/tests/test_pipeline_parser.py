import os
import pytest

import assemblerflow.generator.pipeline_parser as ps
from assemblerflow.tests.data_pipelines import pipelines as pipes


def test_get_lanes():

    raw_string = [
        "A | B)",
        "A | B C D | E F)",
        "A Z | B C (D | E) | G H)",
        "A | B (C | D) | E (E | F I ))"
    ]

    expected = [
        [["A"], ["B"]],
        [["A"], ["B", "C", "D"], ["E", "F"]],
        [["A", "Z"], ["B", "C"], ["G", "H"]],
        [["A"], ["B"], ["E"]]
    ]

    for p, exp in zip(raw_string, expected):
        res = ps.get_lanes(p)
        assert exp == res


def test_linear_connection():

    p = ["A", "B", "C"]
    lane = 1

    res = ps.linear_connection(p, lane)

    assert res == [{
        "input": {
            "process": "A",
            "lane": lane
        },
        "output": {
            "process": "B",
            "lane": lane
        }},
        {"input": {
            "process": "B",
            "lane": lane
        },
        "output": {
            "process": "C",
            "lane": lane
        }
    }]


def test_two_fork_connection():

    source_lane = 1

    res = ps.fork_connection(
        source="A",
        sink=["B", "C"],
        source_lane=source_lane,
        lane=source_lane
    )

    assert res == [{
        "input": {
            "process": "A",
            "lane": source_lane
        },
        "output": {
            "process": "B",
            "lane": source_lane + 1
        }}, {
        "input": {
            "process": "A",
            "lane": source_lane,
        },
        "output": {
            "process": "C",
            "lane": source_lane + 2
        }
    }]


def test_two_fork_connection_mismatch_lane():

    source_lane = 1
    lane = 3

    res = ps.fork_connection(
        source="A",
        sink=["B", "C"],
        source_lane=source_lane,
        lane=lane
    )

    assert res == [{
        "input": {
            "process": "A",
            "lane": source_lane
        },
        "output": {
            "process": "B",
            "lane": lane + 1
        }}, {
        "input": {
            "process": "A",
            "lane": source_lane,
        },
        "output": {
            "process": "C",
            "lane": lane + 2
        }
    }]


def test_multi_fork_connection():

    source_lane = 1

    res = ps.fork_connection(
        source="A",
        sink=["B", "C", "D"],
        source_lane=source_lane,
        lane=source_lane
    )

    assert res == [{
        "input": {
            "process": "A",
            "lane": source_lane
        },
        "output": {
            "process": "B",
            "lane": source_lane + 1
        }}, {
        "input": {
            "process": "A",
            "lane": source_lane,
        },
        "output": {
            "process": "C",
            "lane": source_lane + 2
        }}, {
        "input": {
            "process": "A",
            "lane": source_lane,
        },
        "output": {
            "process": "D",
            "lane": source_lane + 3
        }
    }]


def test_linear_lane_connection():

    res = ps.linear_lane_connection([["A", "B", "C"]], lane=1)

    assert res == [{
        "input": {
            "process": "A",
            "lane": 2
        },
        "output": {
            "process": "B",
            "lane": 2
        }},
        {"input": {
            "process": "B",
            "lane": 2
        },
        "output": {
            "process": "C",
            "lane": 2
        }
    }]


def test_linear_multi_lane_connection():

    res = ps.linear_lane_connection([["A", "B"], ["C", "D"]], lane=1)

    assert res == [{
        "input": {
            "process": "A",
            "lane": 2
        },
        "output": {
            "process": "B",
            "lane": 2
        }},
        {"input": {
            "process": "C",
            "lane": 3
        },
        "output": {
            "process": "D",
            "lane": 3
        }
    }]


def test_get_source_lane():

    pipeline_list = [{'input': {'process': '__init__', 'lane': 1},
                      'output': {'process': 'integrity_coverage', 'lane': 1}},
                     {'input': {'process': 'integrity_coverage', 'lane': 1},
                      'output': {'process': 'fastqc_trimmomatic', 'lane': 1}},
                     {'input': {'process': 'fastqc_trimmomatic', 'lane': 1},
                      'output': {'process': 'spades', 'lane': 2}},
                     {'input': {'process': 'fastqc_trimmomatic', 'lane': 1},
                      'output': {'process': 'skesa', 'lane': 3}}]

    for p, lane in zip(["fastqc_trimmomatic", "spades", "skesa"], [1, 2, 3]):
        res = ps.get_source_lane(p, pipeline_list)
        assert res == lane


def test_parse_pipeline():

    for p, expected in pipes:
        res = ps.parse_pipeline(p)
        assert res == expected


def test_parse_pipeline_file():

    for i in range(1, 9):

        p_path = os.path.join("assemblerflow", "tests", "pipeline_tests",
                              "pipe{}.txt".format(i))
        expected = pipes[i - 1][1]
        print(p_path)
        res = ps.parse_pipeline(p_path)
        assert res == expected
