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

#
plasmid_length = {
    "ACC1": "10",
    "ACC2": "10"
}

# generate a test file
with open("test_depth_file.txt", "w") as text_file:
    text_file.write("seq1\t1\t30")


def test_depth_file_reader():
    """
    test the output of depth_file_reader_function to be a given dict
    """

    # execute the function to get the dict
    result_dict = mapping2json.depth_file_reader(open("test_depth_file.txt"))

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
