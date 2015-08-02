__author__ = 'ddeconti'

import numpy
import pickle
from rdkit.Chem import AllChem, MolFromSmiles
from rdkit import DataStructs
import sys


def parse_smiles(rf, fda_filename):
    try:
        handle = open(fda_filename, 'rU')
    except IOError:
        sys.stderr.write("IOError\n")
        sys.exit()
    fpm_list = []
    line_list = []
    for line in handle:
        arow = line.strip('\n').split('\t')
        smiles = arow[1]
        try:
            mol = MolFromSmiles(smiles)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            arr = numpy.zeros((1,))
            DataStructs.ConvertToNumpyArray(fp, arr)
            fpm = arr
            fpm_list.append(fpm)
        except:
            continue
        if rf.predict(fpm)[0] == 1:
            proba = rf.predict_proba(fpm)[0][1]
            line_list.append(line  + "\t" + str(proba))
    for line in line_list:
        sys.stdout.write(line + "\n")


def main(sa):
    rf_filename = sa[0]
    fda_filename = sa[1]
    rf = pickle.load(open(rf_filename, 'rb'))
    parse_smiles(rf, fda_filename)


if __name__ == "__main__":
    main(sys.argv[1:])
