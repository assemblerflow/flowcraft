import os
import sys
import shutil
import pytest

import flowcraft.flowcraft as af


@pytest.fixture
def tmp():

    os.mkdir("temp")
    yield "temp"
    shutil.rmtree("temp")


def test_check():

    sys.argv.append(1)
    args = af.get_args(["build", "-t 'A B C'", "-c", "-o teste.nf"])

    with pytest.raises(SystemExit):
        af.build(args)


def test_check_invalid():

    sys.argv.append(1)
    args = af.get_args(["build", "-t",  "'A B C()'", "-c", "-o teste.nf"])

    with pytest.raises(SystemExit):
        af.build(args)


def test_build_file(tmp):

    p = os.path.join(os.path.abspath(tmp), "teste.nf")
    sys.argv.append(1)

    args = af.get_args(["build", "-t", "integrity_coverage fastqc", "-o",
                        "{}".format(p)])
    af.build(args)


def test_build_file_2(tmp):

    sys.argv.append(1)
    p = os.path.join(os.path.abspath(tmp), "teste.nf")

    args = af.get_args(["build", "-t integrity_coverage fastqc", "-o",
                        "{}".format(p), "--pipeline-only"])
    af.build(args)

    assert sorted(os.listdir(tmp)) == [".treeDag.json", "containers.config",
                                       "lib", "params.config",
                                       "resources.config", "teste.html",
                                       "teste.nf", "user.config"]
