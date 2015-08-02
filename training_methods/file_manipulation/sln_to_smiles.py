__author__ = 'ddeconti'

import sys
import FileHandler

from rdkit.Chem import MolToSmiles

def print_mol(mol_list):
    for mol in mol_list:
        print MolToSmiles(mol)



def main(sa):
    sln_filename = sa[0]
    sln_list = FileHandler.SlnFile(sln_filename).get_mol_list()
    print_mol(sln_list)


if __name__ == "__main__":
    main(sys.argv[1:])