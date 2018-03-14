import os
import sys
import shutil
import pytest

import assemblerflow.assemblerflow as af


@pytest.fixture
def tmp():

    os.mkdir("temp")
    yield "temp"
    shutil.rmtree("temp")


def test_empty_args():

    args = af.get_args("")

    assert af.check_arguments(args) is False


def test_invalid_path():

    args = af.get_args(["-o invalid/teste.nf"])

    assert af.check_arguments(args) is False


def test_valid_path():

    args = af.get_args(["-o teste.nf"])

    assert af.check_arguments(args) is False


def test_list_short():

    args = af.get_args(["-l"])

    with pytest.raises(SystemExit):
        af.run(args)


def test_list_long():

    args = af.get_args(["-L"])

    with pytest.raises(SystemExit):
        af.run(args)


def test_check():

    args = af.get_args(["-t 'A B C'", "-c", "-o teste.nf"])
    sys.argv.append(1)

    with pytest.raises(SystemExit):
        af.run(args)


def test_build_file(tmp):

    p = os.path.join(os.path.abspath(tmp), "teste.nf")
    sys.argv.append(1)

    args = af.get_args(["-t integrity_coverage fastqc", "-o", "{}".format(p)])
    af.run(args)

    assert os.listdir(tmp) == ["teste.nf"]


def test_build_file(tmp):

    p = os.path.join(os.path.abspath(tmp), "teste.nf")
    sys.argv.append(1)

    args = af.get_args(["-t integrity_coverage fastqc", "-o", "{}".format(p),
                        "--include-templates"])
    af.run(args)

    print(os.listdir(tmp))

    assert os.listdir(tmp) == \
        ["teste.nf", "bin", "templates", "nextflow.config", "lib"]
