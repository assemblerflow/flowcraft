import os
import shutil
import pytest
import flowcraft.templates.integrity_coverage as ic


MAGIC_DICT = {
    b"\x1f\x8b\x08": "gz",
    b"\x42\x5a\x68": "bz2",
    b"\x50\x4b\x03\x04": "zip"
}


def test_no_compression():

    fastq_path = "flowcraft/tests/data/sample.fastq"

    assert ic.guess_file_compression(fastq_path) is None


def test_gzip_compression():

    fastq_path = "flowcraft/tests/data/sample_1.fastq.gz"

    assert ic.guess_file_compression(fastq_path, MAGIC_DICT) == "gz"


def test_bz2_compression():

    fastq_path = "flowcraft/tests/data/sample.fastq.bz2"

    assert ic.guess_file_compression(fastq_path, MAGIC_DICT) == "bz2"


def test_zip_compression():

    fastq_path = "flowcraft/tests/data/sample.fastq.zip"

    assert ic.guess_file_compression(fastq_path, MAGIC_DICT) == "zip"


def test_max_qual_range():

    qual_str = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ' \
               '[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

    assert ic.get_qual_range(qual_str) == (33, 126)


def test_sanger_range():

    sanger_quals = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"

    assert ic.get_qual_range(sanger_quals) == ic.RANGES["Sanger"][1]


def test_solexa_range():

    solexa_quals = ";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"

    assert ic.get_qual_range(solexa_quals) == ic.RANGES["Solexa"][1]


def test_illumina13_range():

    illumina13_quals = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"

    assert ic.get_qual_range(illumina13_quals) == ic.RANGES["Illumina-1.3"][1]


def test_illumina15_range():

    illumina15_quals = "BCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi"

    assert ic.get_qual_range(illumina15_quals) == ic.RANGES["Illumina-1.5"][1]


def test_illumina18_range():

    illumina18_quals = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"

    assert ic.get_qual_range(illumina18_quals) == ic.RANGES["Illumina-1.8"][1]


def test_get_encoding_solexa():

    solexa_quals = "=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"

    quals = ic.get_qual_range(solexa_quals)

    assert ic.get_encodings_in_range(*quals) == (["Solexa"], [64])


def test_get_encoding_illumina18():

    illumina18_quals = "&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"

    quals = ic.get_qual_range(illumina18_quals)

    assert ic.get_encodings_in_range(*quals) == (["Illumina-1.8"], [33])


def test_get_ambiguous_encoding():

    sanger_illumina = "&'()*+,-./0123456789:;<=>?@ABCDEFGHI"

    quals = ic.get_qual_range(sanger_illumina)

    assert ic.get_encodings_in_range(*quals) == (["Sanger", "Illumina-1.8"],
                                                 [33, 33])


@pytest.fixture
def temp_env():
    # Create a temporary directory and move there
    tmp_dir = ".temp"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
    os.chdir(tmp_dir)
    yield
    # Exit and remove temporary directory
    os.chdir(os.pardir)
    shutil.rmtree(tmp_dir)


def test_full_run_files(temp_env):

    opts = [
        "sample",
        ["../flowcraft/tests/data/sample_1.fastq.gz",
         "../flowcraft/tests/data/sample_2.fastq.gz"],
        2.1,
        15,
        ""
    ]

    ic.MAGIC_DICT = MAGIC_DICT
    ic.main(*opts)

    assert sorted(os.listdir(".")) == [
        '.fail', '.report.json', '.status', '.versions', 'sample_coverage',
        'sample_encoding', 'sample_max_len','sample_phred', 'sample_report'
    ]


def test_full_run_success(temp_env):

    opts = [
        "sample",
        ["../flowcraft/tests/data/sample_1.fastq.gz",
         "../flowcraft/tests/data/sample_2.fastq.gz"],
        2.1,
        15,
        ""
    ]

    ic.MAGIC_DICT = MAGIC_DICT
    ic.main(*opts)

    with open("sample_encoding") as fh:
        res = sorted(fh.read().split(","))

    assert res == ["Illumina-1.8", "Sanger"]


def test_low_coverage(temp_env):

    opts = [
        "sample",
        ["../flowcraft/tests/data/sample_1.fastq.gz", "../flowcraft/tests/data/sample_2.fastq.gz"],
        2.1,
        15,
        ""
    ]

    ic.MAGIC_DICT = MAGIC_DICT
    ic.main(*opts)

    with open("sample_coverage") as fh:
        res = fh.read()

    assert res == "fail"


def test_coverage_pass(temp_env):

    opts = [
        "sample",
        ["../flowcraft/tests/data/sample_1.fastq.gz", "../flowcraft/tests/data/sample_2.fastq.gz"],
        2.1e-7,
        15,
        ""
    ]

    ic.MAGIC_DICT = MAGIC_DICT
    ic.main(*opts)

    with open("sample_coverage") as fh:
        res = fh.read()

    assert res == "1438.1"
