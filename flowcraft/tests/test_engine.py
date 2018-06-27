import os
import shutil
import pytest

import flowcraft.generator.engine as eg
import flowcraft.generator.process as pc
import flowcraft.generator.error_handling as eh

from flowcraft.generator.engine import process_map


@pytest.fixture
def single_con():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc", "lane": 1}}
           ]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def single_status():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "skesa", "lane": 1}}]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def single_con_fasta():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "abricate", "lane": 1}}]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def single_con_multi_raw():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "assembly_mapping", "lane": 1}},
           {"input": {"process": "assembly_mapping", "lane": 1},
            "output": {"process": "pilon", "lane": 1}}]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def implicit_link():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc", "lane": 1}},
           {"input": {"process": "fastqc", "lane": 1},
            "output": {"process": "spades", "lane": 1}},
           {"input": {"process": "spades", "lane": 1},
            "output": {"process": "assembly_mapping", "lane": 1}}]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def implicit_link_2():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "spades", "lane": 1}},
           {"input": {"process": "spades", "lane": 1},
            "output": {"process": "assembly_mapping", "lane": 1}}]

    return eg.NextflowGenerator(con, "teste.nf")


@pytest.fixture
def single_fork():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "spades", "lane": 2}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "skesa", "lane": 3}},
           {'input': {'process': 'spades', 'lane': 2},
            'output': {'process': 'abricate', 'lane': 2}},
           {'input': {'process': 'skesa', 'lane': 3},
            'output': {'process': 'abricate', 'lane': 3}}]

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
            "output": {"process": "integrity_coverage", "lane": 3}},
           {"input": {"process": "integrity_coverage", "lane": 3},
            "output": {"process": "check_coverage", "lane": 3}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "spades", "lane": 4}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "skesa", "lane": 5}},
           {"input": {"process": "check_coverage", "lane": 3},
            "output": {"process": "spades", "lane": 6}},
           {"input": {"process": "check_coverage", "lane": 3},
            "output": {"process": "skesa", "lane": 7}}]

    os.mkdir(".temp")
    yield eg.NextflowGenerator(con, os.path.join(".temp", "teste.nf"))
    shutil.rmtree(".temp")


def test_simple_init():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "{}", "lane": 1}}]

    for p in process_map:

        con[0]["output"]["process"] = p
        nf = eg.NextflowGenerator(con, "teste/teste.nf",
                                  ignore_dependencies=True)

        assert [len(nf.processes), nf.processes[1].template] == \
            [2, p]


def test_invalid_process():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "invalid", "lane": 1}}]

    with pytest.raises(SystemExit):
        eg.NextflowGenerator(con, "teste.nf")


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


def test_connections_ignore_type():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "skesa", "lane": 1}},
           {"input": {"process": "skesa", "lane": 1},
            "output": {"process": "patho_typing", "lane": 1}}
           ]

    eg.NextflowGenerator(con, "teste.nf")


def test_build_header(single_con):

    single_con._build_header()

    assert single_con.template != ""


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


def test_connections_channel_update_wfork_2(single_fork):

    p1 = single_fork.processes[3]
    p2 = single_fork.processes[5]

    assert p1.output_channel == p2.input_channel


def test_connections_channel_update_wfork_3(single_fork):

    p1 = single_fork.processes[2]
    p2 = single_fork.processes[4]

    assert p1.output_channel == p2.input_channel



def test_set_channels_single_con_raw_fastq(single_con):

    single_con._set_channels()

    assert [list(single_con.main_raw_inputs.keys())[0],
            len(single_con.main_raw_inputs),
            list(single_con.main_raw_inputs.values())[0]["raw_forks"]] == \
           ["fastq", 1, ["integrity_coverage_in_1_0"]]


def test_set_channels_single_con_raw_fasta(single_con_fasta):

    single_con_fasta._set_channels()

    assert [list(single_con_fasta.main_raw_inputs.keys())[0],
            len(single_con_fasta.main_raw_inputs),
            list(single_con_fasta.main_raw_inputs.values())[0][
                "raw_forks"]] == \
           ["fasta", 1, ["abricate_in_1_0"]]


def test_set_channels_multi_raw_input(single_con_multi_raw):

    single_con_multi_raw._set_channels()

    print(single_con_multi_raw.main_raw_inputs)

    assert [list(single_con_multi_raw.main_raw_inputs.keys()),
            len(single_con_multi_raw.main_raw_inputs)] == \
           [["fasta", "fastq"], 2]


def test_set_channels_secondary_channels_nolink(single_con):

    single_con._set_channels()

    assert single_con.secondary_channels["SIDE_phred"][1]["end"] == []


def test_set_channels_secondary_chanels_link(multi_forks):

    multi_forks._set_channels()

    assert [multi_forks.secondary_channels["SIDE_phred"][1]["end"],
            multi_forks.secondary_channels["SIDE_max_len"][1]["end"],
            multi_forks.secondary_channels["SIDE_max_len"][3]["end"]] == \
           [[], ["SIDE_max_len_4_5"], ["SIDE_max_len_6_7"]]


def test_set_secondary_inputs_raw_forks(raw_forks):

    raw_forks._set_channels()
    raw_forks._set_init_process()

    p = raw_forks

    assert p.main_raw_inputs["fastq"]["raw_forks"] == \
           ["integrity_coverage_in_0_0",
            "patho_typing_in_0_2",
            "seq_typing_in_0_3"]


def test_set_secondary_inputs_multi_raw(single_con_multi_raw):

    single_con_multi_raw._set_channels()
    single_con_multi_raw._set_init_process()

    p = single_con_multi_raw

    assert sorted(list(p.main_raw_inputs.keys())) == ["fasta", "fastq"]


def test_set_secondary_channels(multi_forks):

    multi_forks._set_channels()
    multi_forks._set_secondary_channels()

    p = multi_forks.processes[1]

    print(multi_forks.main_raw_inputs)

    assert [p._context["output_channel"], p._context["forks"]] == \
        ["_integrity_coverage_out_1_0",
         "\n_integrity_coverage_out_1_0.into{ integrity_coverage_out_1_0;"
         "spades_in_1_4;skesa_in_1_5 }\n\n\nSIDE_max_len_1_1.set{"
         " SIDE_max_len_4_5 }\n"]


def test_set_secondary_channels_2(multi_forks):

    multi_forks._set_channels()
    multi_forks._set_secondary_channels()

    p = multi_forks.processes[4]

    assert [p._context["output_channel"], p.main_forks] == \
           ["_check_coverage_out_3_3",
            ["check_coverage_out_3_3", "spades_in_3_6", "skesa_in_3_7"]]


def test_set_implicit_link(implicit_link):

    implicit_link._set_channels()
    implicit_link._set_secondary_channels()

    p = implicit_link.processes[2]

    assert p.main_forks == ["fastqc_out_1_1", "_LAST_fastq_4"]


def test_set_implicit_link(implicit_link_2):

    implicit_link_2._set_channels()
    implicit_link_2._set_secondary_channels()

    p = implicit_link_2.processes[1]

    assert p.main_forks == ["integrity_coverage_out_1_0", "_LAST_fastq_1_3"]


def test_set_status_channels_multi(single_con):

    single_con._set_channels()
    single_con._set_status_channels()

    p = [x for x in single_con.processes[::-1]
         if isinstance(x, pc.StatusCompiler)][0]

    assert p._context["compile_channels"] == \
        "STATUS_integrity_coverage_1_1.mix(STATUS_fastqc2_1_2," \
        "STATUS_fastqc2_report_1_2)"


def test_set_status_channels_single(single_status):

    single_status._set_channels()
    single_status._set_status_channels()

    p = [x for x in single_status.processes[::-1]
         if isinstance(x, pc.StatusCompiler)][0]

    assert p._context["compile_channels"] == "STATUS_skesa_1_1"


def test_set_compiler_channels(single_status):

    single_status.lane = 1
    single_status._set_channels()
    single_status._set_compiler_channels()

    p = [x for x in single_status.processes[::-1]
         if isinstance(x, pc.StatusCompiler)][0]

    assert p._context["compile_channels"] == "STATUS_skesa_1_1"


def test_set_status_channels_no_status(single_status):

    single_status.processes[1].status_channels = []

    single_status._set_channels()
    single_status._set_status_channels()

    with pytest.raises(IndexError):
        p = [x for x in single_status.processes[::-1]
             if isinstance(x, pc.StatusCompiler)][0]


def test_set_status_channels_duplicate_status(single_status):

    single_status.processes[1].status_channels = ["A", "A"]

    single_status._set_channels()

    with pytest.raises(eh.ProcessError):
        single_status._set_status_channels()


def test_build(multi_forks):

    multi_forks.build()

    assert multi_forks.template != ""


def test_resources_string(single_con):

    res_dict = {"procA": {"cpus": 1, "memory": "'4GB'", "container": "img",
                          "version": "1"}}

    res = single_con._get_resources_string(res_dict, 1)

    assert res == '\n\t$procA_1.cpus = 1\n\t$procA_1.memory = \'4GB\''


def test_resources_string_2(single_con):

    res_dict = {"procA": {"cpus": 1, "container": "img",
                          "version": "1"}}

    res = single_con._get_resources_string(res_dict, 1)

    assert res == '\n\t$procA_1.cpus = 1'


def test_resources_string_3(single_con):

    res_dict = {"procA": {"cpus": 1, "memory": "'4GB'", "container": "img",
                          "version": "1"},
                "procB": {"memory": "{ 4.GB * task.attempt }"}}

    res = single_con._get_resources_string(res_dict, 1)

    assert res == '\n\t$procA_1.cpus = 1\n\t$procA_1.memory = \'4GB\'' \
                  '\n\t$procB_1.memory = { 4.GB * task.attempt }'


def test_container_string(single_con):

    res_dict = {"procA": {"cpus": 1, "memory": "4GB", "container": "img",
                          "version": "1"}}

    res = single_con._get_container_string(res_dict, 2)

    assert res == '\n\t$procA_2.container = "img:1"'


def test_container_string_2(single_con):

    res_dict = {"procA": {"cpus": 1, "memory": "4GB", "container": "img",
                          "version": "1"},
                "procB": {"container": "img"}}

    res = single_con._get_container_string(res_dict, 2)

    assert res == '\n\t$procA_2.container = "img:1"\n\t' \
                  '$procB_2.container = "img:latest"'


def test_extra_inputs_1():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc={'extra_input':'teste'}", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert nf.processes[2].extra_input == "teste"


def test_extra_inputs_2():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "spades", "lane": 1}},
           {"input": {"process": "spades", "lane": 1},
            "output": {"process": "abricate={'extra_input':'teste'}", "lane": 1}}
           ]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert nf.processes[3].extra_input == "teste"


def test_extra_inputs_3():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc={'extra_input':'teste'}", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")
    nf._set_channels()

    assert [list(nf.extra_inputs.keys())[0],
            nf.extra_inputs["teste"]["input_type"],
            nf.extra_inputs["teste"]["channels"]] == \
           ["teste", "fastq", ["EXTRA_fastqc_1_2"]]


def test_extra_inputs_default():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "spades", "lane": 1}},
           {"input": {"process": "spades", "lane": 1},
            "output": {"process": "abricate={'extra_input':'default'}", "lane": 1}}
           ]

    nf = eg.NextflowGenerator(con, "teste.nf")
    nf._set_channels()

    assert [list(nf.extra_inputs.keys())[0],
            nf.extra_inputs["fasta"]["input_type"],
            nf.extra_inputs["fasta"]["channels"]] == \
           ["fasta", "fasta", ["EXTRA_abricate_1_3"]]


def test_extra_inputs_invalid():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc={'extra_input':'default'}", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    with pytest.raises(SystemExit):
        nf._set_channels()


def test_extra_inputs_invalid_2():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "spades={'extra_input':'teste'}", "lane": 1}},
           {"input": {"process": "spades", "lane": 1},
            "output": {"process": "abricate={'extra_input':'teste'}", "lane": 1}}
           ]

    nf = eg.NextflowGenerator(con, "teste.nf")

    with pytest.raises(SystemExit):
        nf._set_channels()


def test_run_time_directives():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc={'cpus':'3'}", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert nf.processes[2].directives["fastqc2"]["cpus"] == "3"


def test_run_time_directives_full():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc={'cpus':'3','memory':'4GB',"
                                  "'container':'img','version':'1'}",
                       "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert [nf.processes[2].directives["fastqc2"]["cpus"],
            nf.processes[2].directives["fastqc2"]["memory"],
            nf.processes[2].directives["fastqc2"]["container"],
            nf.processes[2].directives["fastqc2"]["version"]] == \
           ["3", "4GB", "img", "1"]


def test_run_time_directives_invalid():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "integrity_coverage", "lane": 1}},
           {"input": {"process": "integrity_coverage", "lane": 1},
            "output": {"process": "fastqc={'cpus'", "lane": 1}}]

    with pytest.raises(SystemExit):
        eg.NextflowGenerator(con, "teste.nf")


def test_not_automatic_dependency():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "spades", "lane": 1}}]

    with pytest.raises(SystemExit):
        eg.NextflowGenerator(con, "teste.nf", auto_dependency=False)


def test_automatic_dependency():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "spades", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert nf.processes[1].template == "integrity_coverage"


def test_automatic_dependency_2():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "spades", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert nf.processes[1].output_channel == nf.processes[2].input_channel


def test_automatic_dependency_3():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "spades", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert [nf.processes[1].parent_lane, nf.processes[2].parent_lane] == \
           [None, 1]


def test_automatic_dependency_wfork():

    con = [{"input": {"process": "__init__", "lane": 0},
            "output": {"process": "spades", "lane": 1}},
           {"input": {"process": "__init__", "lane": 0},
            "output": {"process": "integrity_coverage", "lane": 2}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert nf.processes[1].template == "integrity_coverage"


def test_automatic_dependency_wfork_2():

    con = [{"input": {"process": "__init__", "lane": 0},
            "output": {"process": "spades", "lane": 1}},
           {"input": {"process": "__init__", "lane": 0},
            "output": {"process": "integrity_coverage", "lane": 2}}]

    nf = eg.NextflowGenerator(con, "teste.nf")
    nf._set_channels()

    assert len(nf.main_raw_inputs["fastq"]["raw_forks"]) == 2


def test_automatic_dependency_wfork_3():

    con = [{"input": {"process": "__init__", "lane": 0},
            "output": {"process": "reads_download", "lane": 1}},
           {"input": {"process": "reads_download", "lane": 1},
            "output": {"process": "skesa", "lane": 2}},
           {"input": {"process": "reads_download", "lane": 1},
            "output": {"process": "spades", "lane": 3}}
           ]

    nf = eg.NextflowGenerator(con, "teste.nf")
    nf._set_channels()

    assert nf.processes[3].parent_lane == 1


def test_automatic_dependency_wfork_4():

    con = [{"input": {"process": "__init__", "lane": 0},
            "output": {"process": "reads_download", "lane": 1}},
           {"input": {"process": "reads_download", "lane": 1},
            "output": {"process": "skesa", "lane": 2}},
           {"input": {"process": "reads_download", "lane": 1},
            "output": {"process": "spades", "lane": 3}}
           ]

    nf = eg.NextflowGenerator(con, "teste.nf")
    nf._set_channels()

    assert nf.processes[4].parent_lane == 3


def test_automatic_dependency_multi():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "trimmomatic", "lane": 1}},
           {"input": {"process": "trimmomatic", "lane": 1},
            "output": {"process": "spades", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert len([x for x in nf.processes
                if x.template == "integrity_coverage"]) == 1


def test_automatic_dependency_non_raw():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "spades", "lane": 1}},
           {"input": {"process": "spades", "lane": 1},
            "output": {"process": "pilon", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    assert nf.processes[2].parent_lane == 1


def test_patlas_compiler_channels():

    con = [{"input": {"process": "__init__", "lane": 0},
            "output": {"process": "mash_screen", "lane": 1}},
           {"input": {"process": "__init__", "lane": 0},
            "output": {"process": "mapping_patlas", "lane": 2}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    nf._set_channels()
    nf._set_compiler_channels()

    assert len(nf.compilers["patlas_consensus"]["channels"]) == 2


def test_patlas_compiler_channels_2():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "mash_screen", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    nf._set_channels()
    nf._set_compiler_channels()

    assert len(nf.compilers["patlas_consensus"]["channels"]) == 1


def test_patlas_compiler_channels_empty():

    con = [{"input": {"process": "__init__", "lane": 1},
            "output": {"process": "trimmomatic", "lane": 1}}]

    nf = eg.NextflowGenerator(con, "teste.nf")

    nf._set_channels()
    nf._set_compiler_channels()

    assert len(nf.compilers["patlas_consensus"]["channels"]) == 0
