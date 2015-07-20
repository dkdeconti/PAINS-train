__author__ = 'ddeconti'

import sys

def parse_chemicals(filename):
    chem_dict = {}
    try:
        handle = open(filename, 'rU')
    except IOError as e:
        sys.stderr.write("IOError: " + str(e) +
                         "\nError in parse_chemicals()\n")
        sys.exit()
    for line in handle:
        line = line.strip('\n').split('\t')
        chem_id = line[0]
        smile = line[-1]
        chem_dict[chem_id] = {"smile": smile, "num_inter": 0}
    handle.close()
    return chem_dict

def parse_interations(filename, chem_dict):
    try:
        handle = open(filename, 'rU')
    except IOError as e:
        sys.stderr.write("IOError: " + str(e) +
                         "\nError in parse_interactions()\n")
        sys.exit()
    for line in handle:
        line = line.strip('\n').split('\t')
        chem_id = line[0]
        if chem_id in chem_dict:
            chem_dict[chem_id]["num_inter"] += 1
    handle.close()
    return chem_dict


def main(sa):
    chems_filename = sa[0]
    interactions_filename = sa[1]


if __name__ == "__main__":
    main(sys.argv[1:])