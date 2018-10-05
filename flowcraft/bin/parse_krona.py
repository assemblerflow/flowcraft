#!/usr/bin/env python2

import json
import csv
import sys



def main(file_path):


    # Parse krona TSV results file
    file_parsed = list(csv.reader(open(file_path, 'r'), delimiter='\t'))

    dict_keys = {}
    """
    dict: Used to count the total number of results for each key level
    """
    levels = []
    """
    list: List of the resulting levels, with values, keys, and children
    """

    for line in file_parsed:
        value = 0
        # Count total number of results for each level. If line only has a
        # value, it is classifier as Unclassified
        for index, level in enumerate(line):
            if len(line) == 1:
                level_to_use = "Unclassified"
                value = level
            else:
                level_to_use = level

            if index > 0 or len(line) == 1:
                if level_to_use not in dict_keys:
                    dict_keys[level_to_use] = float(value)
                else:
                    dict_keys[level_to_use] += float(value)
            else:
                value = level_to_use

    for line in file_parsed:
        current_level = {}
        # Reads the array reversed and creates the dictionary tree of
        # classification of children. Last element has and empty object has
        # children.
        for index, level in enumerate(line[::-1]):
            try:
                current_level = {
                    "key": level,
                    "value": dict_keys[level],
                    "children": current_level
                }
            except KeyError:
                print("Value not in dictionary")

        levels.append(current_level)

    with open("parsed_krona.txt", "w") as k:
        k.write(json.dumps(levels))


if __name__ == "__main__":
    main(sys.argv[1])

