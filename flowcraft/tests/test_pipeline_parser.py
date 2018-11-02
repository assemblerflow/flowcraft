import os
import json

import flowcraft.generator.pipeline_parser as ps
from flowcraft.tests.data_pipelines import pipelines as pipes


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

    res = ps.get_source_lane(["integrity_coverage", "fastqc_trimmomatic"],
                             pipeline_list)

    assert res == 1


def test_get_source_lane_2():

    pipeline_list = [{'input': {'process': '__init__', 'lane': 1},
                      'output': {'process': 'integrity_coverage', 'lane': 1}},
                     {'input': {'process': 'integrity_coverage', 'lane': 1},
                      'output': {'process': 'fastqc_trimmomatic', 'lane': 1}},
                     {'input': {'process': 'fastqc_trimmomatic', 'lane': 1},
                      'output': {'process': 'spades', 'lane': 2}},
                     {'input': {'process': 'fastqc_trimmomatic', 'lane': 1},
                      'output': {'process': 'skesa', 'lane': 3}},
                     {'input': {'process': 'spades', 'lane': 2},
                      'output': {'process': 'pilon', 'lane': 2}},
                     {'input': {'process': 'skesa', 'lane': 3},
                      'output': {'process': 'pilon', 'lane': 3}},
                     ]

    res = ps.get_source_lane(["spades", "pilon"], pipeline_list)

    assert res == 2


def test_parse_pipeline():

    for p, expected in pipes:
        res = ps.parse_pipeline(p)
        assert res == expected


def test_parse_pipeline_file():

    for i in range(1, 9):

        p_path = os.path.join("flowcraft", "tests", "pipeline_tests",
                              "pipe{}.txt".format(i))
        expected = pipes[i - 1][1]
        print(p_path)
        res = ps.parse_pipeline(p_path)
        print(res)
        assert res == expected


def test_unique_id_len():

    pip_list = [
        "A B C",
        "A (B C (D | E)| B C (D | E))",
        "A (B C (D | E)| C (D | E))",
        "A (B C (D | E)| B (D | E))",
    ]

    res_list = [
        "A_0 B_1 C_2",
        "A_0 (B_1 C_2 (D_3 | E_4)| B_5 C_6 (D_7 | E_8))",
        "A_0 (B_1 C_2 (D_3 | E_4)| C_5 (D_6 | E_7))",
        "A_0 (B_1 C_2 (D_3 | E_4)| B_5 (D_6 | E_7))",
    ]

    for x, pip_str in enumerate(pip_list):
        res_str, res_ids = ps.add_unique_identifiers(pip_str)
        assert res_str.replace(" ", "") == res_list[x].replace(" ", "")

def test_remove_id():

    pip_list = [
        "A B C",
        "A (B C (D | E)| B C (D | E))",
    ]

    pipeline_mod_links = [
        [{'input': {'process': '__init__', 'lane': 1},
          'output': {'process': 'A_0', 'lane': 1}},
         {'input': {'process': 'A_0', 'lane': 1},
          'output': {'process': 'B_1', 'lane': 1}},
         {'input': {'process': 'B_1', 'lane': 1},
          'output': {'process': 'C_2', 'lane': 1}}],
        [{'input': {'process': '__init__', 'lane': 1},
          'output': {'process': 'A_0', 'lane': 1}},
         {'input': {'process': 'A_0', 'lane': 1},
          'output': {'process': 'B_1', 'lane': 2}},
         {'input': {'process': 'A_0', 'lane': 1},
          'output': {'process': 'B_5', 'lane': 3}},
         {'input': {'process': 'B_1', 'lane': 2},
          'output': {'process': 'C_2', 'lane': 2}},
         {'input': {'process': 'B_5', 'lane': 3},
          'output': {'process': 'C_6', 'lane': 3}},
         {'input': {'process': 'C_2', 'lane': 2},
          'output': {'process': 'D_3', 'lane': 4}},
         {'input': {'process': 'C_2', 'lane': 2},
          'output': {'process': 'E_4', 'lane': 5}},
         {'input': {'process': 'C_6', 'lane': 3},
          'output': {'process': 'D_7', 'lane': 6}},
         {'input': {'process': 'C_6', 'lane': 3},
          'output': {'process': 'E_8', 'lane': 7}}]
    ]

    pipeline_exp_links = [
        [{'input': {'process': '__init__', 'lane': 1},
          'output': {'process': 'A', 'lane': 1}},
         {'input': {'process': 'A', 'lane': 1},
          'output': {'process': 'B', 'lane': 1}},
         {'input': {'process': 'B', 'lane': 1},
          'output': {'process': 'C', 'lane': 1}}],
        [{'input': {'process': '__init__', 'lane': 1},
          'output': {'process': 'A', 'lane': 1}},
         {'input': {'process': 'A', 'lane': 1},
          'output': {'process': 'B', 'lane': 2}},
         {'input': {'process': 'A', 'lane': 1},
          'output': {'process': 'B', 'lane': 3}},
         {'input': {'process': 'B', 'lane': 2},
          'output': {'process': 'C', 'lane': 2}},
         {'input': {'process': 'B', 'lane': 3},
          'output': {'process': 'C', 'lane': 3}},
         {'input': {'process': 'C', 'lane': 2},
          'output': {'process': 'D', 'lane': 4}},
         {'input': {'process': 'C', 'lane': 2},
          'output': {'process': 'E', 'lane': 5}},
         {'input': {'process': 'C', 'lane': 3},
          'output': {'process': 'D', 'lane': 6}},
         {'input': {'process': 'C', 'lane': 3},
          'output': {'process': 'E', 'lane': 7}}]
    ]

    for x, pip_str in enumerate(pip_list):
        res_str, res_ids = ps.add_unique_identifiers(pip_str)
        res = ps.remove_unique_identifiers(res_ids, pipeline_mod_links[x])
        assert json.dumps(res) == json.dumps(pipeline_exp_links[x])