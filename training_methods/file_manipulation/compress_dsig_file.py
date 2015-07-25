__author__ = 'ddeconti'

'''
Compresses total drugs by number of interactions.
'''

import re
import sys


def print_drug(name, num):
    print "Test", name, num
    outstr = '\t'.join([name, str(num)])
    sys.stdout.write(outstr + '\n')


def parse_chem(filename):
    try:
        handle = open(filename, 'rU')
    except IOError as e:
        sys.stderr.write("IOError: " + str(e) +
                         "\nError in parse_chem()\n")
        sys.exit()
    current_drug = None
    current_num = 0
    handle.readline()
    for line in handle:
        line = line.strip('\n').split('\t')
        if not current_drug:
            if re.search("D1", line[3]):
                current_drug = line[0]
                current_num += 1
                continue
            else:
                continue
        if current_drug != line[0]:
            print_drug(current_drug, current_num)
            if re.search("D1", line[3]):
                current_drug = line[0]
                current_num = 1
            else:
                print "Is this the culprit?"
                current_drug = None
                current_num = 0
    if current_drug:
        print "last one:", current_drug, current_num
        print_drug(current_drug, current_num)
    handle.close()


def main(sa):
    chem_filename = sa[0]
    parse_chem(chem_filename)
    pass


if __name__ == "__main__":
    main(sys.argv[1:])