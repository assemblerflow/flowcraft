import os
import pytest

import assemblerflow.generator.process as pc
import assemblerflow.generator.error_handling as eh

from assemblerflow.generator.engine import process_map


@pytest.fixture
def mock_process():

    return pc.Process(template="integrity_coverage")


@pytest.fixture
def process_wchannels():

    p = pc.Process(template="integrity_coverage")

    p.input_channel = "in_channel"
    p.output_channel = "out_channel"

    return p


@pytest.fixture
def mock_status():

    return pc.Status(template="status_compiler")


@pytest.fixture
def mock_init():

    return pc.Init(template="init")


def test_process_init():

    for template, proc in process_map.items():

        p = proc(template=template)

        assert p.template == template


def test_set_correct_template(mock_process):

    mock_process._set_template("fastqc")

    assert os.path.exists(mock_process._template_path)


def test_set_wrong_template(mock_process):

    with pytest.raises(eh.ProcessError):
        mock_process._set_template("wrong_template")


def test_template_render(process_wchannels):

    process_wchannels.set_channels(pid=1)
    t = process_wchannels.template_str

    assert 1


def test_main_channel_setup(mock_process):

    mock_process.set_main_channel_names("input_suf", "output_suf", 1)

    assert [mock_process.input_channel.endswith("input_suf"),
            mock_process.output_channel.endswith("output_suf"),
            mock_process.lane] == [True, True, 1]


def test_main_raw_channel_self(mock_process):
    """Tests the retrieval of the raw input channel when the input type is
    inferred from the class"""

    mock_process.input_type = "fastq"
    res = mock_process.get_user_channel("myChannel")

    assert res == {"input_channel": "myChannel",
                   **mock_process.RAW_MAPPING["fastq"]}


def test_main_raw_channel_fastq(mock_process):

    res = mock_process.get_user_channel("myChannel", "fastq")

    assert res == {"input_channel": "myChannel",
                   **mock_process.RAW_MAPPING["fastq"]}


def test_main_raw_channel_fasta(mock_process):

    res = mock_process.get_user_channel("myChannel", "fasta")

    assert res == {"input_channel": "myChannel",
                   **mock_process.RAW_MAPPING["fasta"]}


def test_main_raw_channel_invalid(mock_process):

    res = mock_process.get_user_channel("myChannel", "invalid")

    assert res is None


def test_channels_setup(process_wchannels):

    process_wchannels.set_channels(pid=1)

    expected = {"input_channel": "in_channel",
                "output_channel": "out_channel",
                "template": process_wchannels.template,
                "pid": 1,
                "forks": ""}

    assert process_wchannels._context == expected


def test_channels_setup_withforks(process_wchannels):

    process_wchannels.forks = ["A", "B"]

    process_wchannels.set_channels(pid=1)

    expected = {"input_channel": "in_channel",
                "output_channel": "out_channel",
                "template": process_wchannels.template,
                "pid": 1,
                "forks": "A\nB"}

    assert process_wchannels._context == expected


def test_setup_one_raw_fork(process_wchannels):

    process_wchannels.main_forks = ["A"]
    process_wchannels.set_channels(pid=1)

    expected = {"input_channel": "in_channel",
                "output_channel": "out_channel",
                "template": process_wchannels.template,
                "pid": 1,
                "forks": "\nout_channel.set{ A }\n"}

    assert process_wchannels._context == expected


def test_setup_multiple_raw_forks(process_wchannels):

    process_wchannels.main_forks = ["A", "B"]
    process_wchannels.set_channels(pid=1)

    expected = {"input_channel": "in_channel",
                "output_channel": "out_channel",
                "template": process_wchannels.template,
                "pid": 1,
                "forks": "\nout_channel.into{ A;B }\n"}

    assert process_wchannels._context == expected


def test_channels_setup_status(process_wchannels):

    process_wchannels.status_channels = ["A", "B"]

    process_wchannels.set_channels(pid=1)

    assert process_wchannels.status_strs == ["A_1", "B_1"]


def test_update_main_fork_noprevious(process_wchannels):
    """Updates the forks attributes when there are no previous main forks"""

    process_wchannels.set_channels(pid=1)
    process_wchannels.update_main_forks("A")

    assert [process_wchannels.output_channel,
            process_wchannels.main_forks,
            process_wchannels.forks] == \
           ["_out_channel",
            ["out_channel", "A"],
            ["\n_out_channel.into{ out_channel;A }\n"]]


def test_secondary_channels_multisink(process_wchannels):

    process_wchannels.set_channels(pid=1)
    process_wchannels.set_secondary_channel("A", ["B", "C"])

    assert process_wchannels.forks == ["\nA_1.into{ B;C }\n"]


def test_secondary_channels_singlesink(process_wchannels):

    process_wchannels.set_channels(pid=1)
    process_wchannels.set_secondary_channel("A", ["B"])

    assert process_wchannels.forks == ["\nA_1.set{ B }\n"]


def test_secondary_channels_duplicatesink(process_wchannels):

    process_wchannels.set_channels(pid=1)
    process_wchannels.set_secondary_channel("A", ["B", "B"])

    assert process_wchannels.forks == ["\nA_1.set{ B }\n"]


def test_status_init(mock_status):

    assert mock_status.template == "status_compiler"


def test_status_channel_setup_empty(mock_status):

    with pytest.raises(eh.ProcessError):
        mock_status.set_status_channels([])


def test_status_channel_single(mock_status):

    mock_status.set_status_channels(["A"])

    assert mock_status._context == {"status_channels": "A"}


def test_status_channel_two(mock_status):

    mock_status.set_status_channels(["A", "B"])

    assert mock_status._context == {"status_channels": "A.mix(B)"}


def test_status_channel_multiple(mock_status):

    mock_status.set_status_channels(["A", "B", "C"])

    assert mock_status._context == {"status_channels": "A.mix(B,C)"}


def test_init_process(mock_init):

    assert mock_init.template == "init"


def test_init_raw_inputs_single(mock_init):

    mock_init.set_raw_inputs({"A": {"channel": "rawChannel",
                                    "raw_forks": ["A"],
                                    "channel_str": "rawChannel.Channel"}})

    assert [mock_init.forks, mock_init._context["main_inputs"]] == \
        [["\nrawChannel.set{ A }\n"], "rawChannel.Channel"]


def test_init_raw_inputs_multi_forks(mock_init):

    mock_init.set_raw_inputs({"A": {"channel": "rawChannel",
                                    "raw_forks": ["A", "B"],
                                    "channel_str": "rawChannel.Channel"}})

    assert [mock_init.forks, mock_init._context["main_inputs"]] == \
        [["\nrawChannel.into{ A;B }\n"], "rawChannel.Channel"]


def test_init_multi_raw_inputs(mock_init):

    mock_init.set_raw_inputs({"A": {"channel": "rawChannel",
                                    "raw_forks": ["A", "B"],
                                    "channel_str": "rawChannel.Channel"},
                              "B": {"channel": "otherChannel",
                                    "raw_forks": ["C"],
                                    "channel_str": "otherChannel.Channel"}})

    assert [mock_init.forks, mock_init._context["main_inputs"]] == \
        [["\nrawChannel.into{ A;B }\n",
          "\notherChannel.set{ C }\n"],
         "rawChannel.Channel\notherChannel.Channel"]


def test_init_secondary_inputs(mock_init):

    mock_init.set_secondary_inputs(
        {"genomeSize": "IN_genome_size = Channel.value(params.genomeSize)"})

    assert mock_init._context["secondary_inputs"] == \
        "IN_genome_size = Channel.value(params.genomeSize)"


def test_init_multi_secondary_inputs(mock_init):

    mock_init.set_secondary_inputs(
        {"genomeSize": "IN_genome_size = Channel.value(params.genomeSize)",
         "other": "Other"})

    assert mock_init._context["secondary_inputs"] == \
        "IN_genome_size = Channel.value(params.genomeSize)\nOther"
