import os
import json
import pytest
import flowcraft.templates.pATLAS_consensus_json as pATLAS_consensus_json


@pytest.fixture
def fetch_file(tmpdir, request):
    # create a temporary file mash screen txt file
    mash_file = tmpdir.join("mash_file.txt")
    depth_file = tmpdir.join("test_depth_file.txt")
    mash_file.write('{"ACC1": ["0.9", "1"]}')
    depth_file.write('{"ACC1": 0.9}')

    list_of_files = [str(mash_file), str(depth_file)]

    pATLAS_consensus_json.main(list_of_files)
    result_dict = json.load(open("consensus_file.json"))

    # finalizer statement that removes .report.json
    def remove_test_files():
        os.remove(".report.json")
    request.addfinalizer(remove_test_files)

    return result_dict


def test_generate_file(fetch_file):
    """
    Check if consensus output file is generated
    """
    assert os.path.isfile("consensus_file.json")


def test_generate_report(fetch_file):
    """
    This tests if the report.json file is generated
    """
    assert os.path.isfile(".report.json")


def test_generated_dict(fetch_file):
    """
    This function checks if the file contains a dict
    """
    result_dict = fetch_file
    assert isinstance(result_dict, dict)


def test_generated_dict_contents1(fetch_file):
    """
    Checks if accession in both expected and result dict is the same
    """
    result_dict = fetch_file
    assert list(result_dict.keys())[0] == "ACC1"


def test_generated_dict_contents2(fetch_file):
    """
    checks if the resulting values for each type of file are the proper type
    """
    result_dict = fetch_file
    first_file_values = list(result_dict[
                                 list(result_dict.keys())[0]].values())[0]
    second_file_values = list(result_dict[
                                  list(result_dict.keys())[0]].values())[1]
    list_of_checks = [isinstance(first_file_values, list),
                      isinstance(second_file_values, float)]

    assert all(list_of_checks)


def test_generated_dict_contents3(fetch_file):
    """
    checks that the accession in dict must have two keys
    """
    result_dict = fetch_file
    assert len(list(result_dict[list(result_dict.keys())[0]].keys())) == 2
