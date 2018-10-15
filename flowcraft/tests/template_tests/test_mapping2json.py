import flowcraft.templates.mapping2json as templates_fc

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


def test_depth_file_reader():
    perc_bases_cov, dict_cov = templates_fc.generate_jsons(depth_dict_coverage,
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

