import os
import json
import pytest
import flowcraft.templates.mapping2json as mapping2json

depth_dict_coverage = {
    "ACC1": {
        "1": 20,
        "2": 30,
        "3": 50,
        "4": 40,
        "5": 30,
        "6": 20,
        "7": 20,
        "8": 40,
        "9": 30,
        "10": 20
    },
    "ACC2": {
        "1": 20,
        "2": 30,
        "3": 50,
        "4": 40,
        "7": 20,
        "8": 40,
        "9": 30,
        "10": 20
    }
}

plasmid_length = {
    "ACC1": "10",
    "ACC2": "10"
}


@pytest.fixture
def fetch_file(tmpdir, request):
    # create a temporary file depth txt file
    depth_file = tmpdir.join("test_depth_file.txt")
    depth_file.write("ACC1\t1\t30")
    # creates a temporary file with the json_dict
    json_dict = tmpdir.join("test_json_dict.json")
    json_dict.write('{"ACC1": 2000}')

    # finalizer statement that removes .report.json
    def remove_test_files():
        os.remove(".report.json")
    request.addfinalizer(remove_test_files)

    return json_dict, str(depth_file)


def test_depth_file_reader(tmpdir):
    """
    test the output of depth_file_reader_function to be a given dict
    """

    # create a temporary file to make a dict from
    file_handle = tmpdir.join("test_depth_file.txt")
    file_handle.write("seq1\t1\t30")

    # execute the function to get the dict
    result_dict = mapping2json.depth_file_reader(open(str(file_handle)))

    assert result_dict == {"seq1": {"1": 30}}


def test_generate_jsons():
    """
    Test if the returns from this function are the expected results
    """
    perc_bases_cov, dict_cov = mapping2json.generate_jsons(depth_dict_coverage,
                                                           plasmid_length, 0.9)

    # asserts if all the returned values are the expected ones
    assert (perc_bases_cov, dict_cov) == (
        {"ACC1": 1.0},
        {"ACC1": {
            "length": 10,
            "interval": 1,
            "values": [20, 30, 50, 40, 30, 20, 20, 40, 30, 20]
        }}
    )


def test_generate_file(fetch_file):
    """
    This function tests if the output json file is generated
    """
    json_dict, depth_file = fetch_file
    # executes the function to be tested
    mapping2json.main(str(depth_file), str(json_dict), "0", "test")
    assert os.path.isfile("{}_mapping.json".format(depth_file))


def test_generate_report(fetch_file):
    """
    This tests if the report.json file is generated
    """
    json_dict, depth_file = fetch_file
    # executes the function to be tested
    mapping2json.main(str(depth_file), str(json_dict), "0", "test")
    assert os.path.isfile(".report.json")


def test_generated_dict(fetch_file):
    """
    This function checks if the file contains a dict
    """
    json_dict, depth_file = fetch_file
    # executes the function to be tested
    mapping2json.main(str(depth_file), str(json_dict), "0", "test")
    result_dict = json.load(open("{}_mapping.json".format(str(depth_file))))
    assert isinstance(result_dict, dict)


def test_generated_dict_contents(fetch_file):
    """
    This function tests if the dictionary in the file has the required fields
    """
    expected_dict = {"ACC1": 0.0}
    json_dict, depth_file = fetch_file
    # executes the function to be tested
    mapping2json.main(str(depth_file), str(json_dict), "0", "test")
    result_dict = json.load(open("{}_mapping.json".format(str(depth_file))))
    assert result_dict == expected_dict
