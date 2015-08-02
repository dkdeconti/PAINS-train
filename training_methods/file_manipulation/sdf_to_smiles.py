__author__ = 'ddeconti'

import sys

from rdkit.Chem import SDMolSupplier, MolToSmiles

def sdf_to_smile(sdf_filename):
    sdf_struct = SDMolSupplier(sdf_filename)
    for mol in sdf_struct:
        try:
            print MolToSmiles(mol)
        except:
            continue

def main(sa):
    sdf_filename = sa[0]
    sdf_to_smile(sdf_filename)
    pass

if __name__ == "__main__":
    main(sys.argv[1:])