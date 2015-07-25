__author__ = 'ddeconti'

import FileObjects
import re
import sys

from rdkit.Chem import AllChem, MolFromSmiles
from rdkit.Chem.rdSLNParse import MolFromSLN
from rdkit.Chem import SDMolSupplier


class WikiScrapedDB():
    '''
    Custom wrapper around drug name - side effect counts - gene count - SMILES
    Produce a dictionary of the values
    '''
    def __init__(self, filename):
        '''
        :param fileanme:
        :return: None

        Creates a container for dict with key values
        '''
        try:
            handle = open(filename, 'rU')
        except IOError as e:
            sys.stderr.write("IOError: " + str(e) +
                             "\nError in WikiScrapedDB.__init__()\n")
            sys.exit()
        self.drug_list = []
        handle.readline()
        for line in handle:
            line = line.strip('\n').split('\t')
            drug_name = line[0]
            try:
                se_count = int(line[1])
                gene_count = int(line[2])
            except ValueError as e:
                sys.stderr.write("ValueError: " + str(e) +
                                 "\nError in WikiSrapedDB.__init__()\n")
                sys.exit()
            smiles = line[3]
            self.drug_list.append(FileObjects.DrugEffects(drug_name, se_count,
                                                          gene_count, smiles))
        handle.close()

    def get_drug_list(self):
        return self.drug_list


class SdfFile():
    '''
    Custom wrapper around RDKit's Struct file input methods.
    Focus is to produce a list of RDKit.Mol object.
    Acts as a container for the RDKit.Mol object list.
    '''
    def __init__(self, filename):
        '''
        :param filename:
        :return: None

        Creates container for reference for RDKit.Mol object list:
            self.mol_list
        '''
        self.filename = filename
        try:
            sdl_struct = SDMolSupplier(self.filename)
        except IOError as e:
            sys.stderr.write("IOError: " + str(e) +
                             "\nError in SdfFile.__init__()\n")
            sys.stderr.flush()
            sys.exit()
        try:
            self.mol_list = [m for m in sdl_struct]
        except:
            sys.stderr.write("Error parsing SDMolSupplier object\n" +
                             "Error in SdfFile.__init__()\n")
            sys.exit()
        self.fingerprint_list = None

    def mol_to_fp(self):
        fp_list = []
        for mol in self.get_mol_list():
            try:
                fp_list.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2))
                #fp_list.append(FingerprintMols.FingerprintMol(mol))
            except:
                continue
        #return filter(lambda x: x != None, fp_list)
        return fp_list

    # Setters

    def set_fingerprint_list(self, fp_list):
        self.fingerprint_list = fp_list

    # Getters

    def get_fingerprint_list(self):
        if not self.fingerprint_list:
            self.set_fingerprint_list(self.mol_to_fp())
        return self.fingerprint_list

    def get_mol_list(self):
        return self.mol_list


class SlnFile():
    '''
    Custom wrapper around RDKit's SLN parsers.
    Focus is to produce a list of RDKit.Mol objects.
    Acts as container for the RDKit.Mol object list.
    '''
    def __init__(self, filename):
        self.filename = filename
        self.sln_list = self.parse_sln_file()
        self.mol_list = None
        self.fingerprint_list = None

    def parse_sln_file(self):
        '''
        Specifially parses sybyl notation (SLN file) from
        JA Holloway 2010 PAINs paper

        :return: list(rdkit.Chem Mol)
        '''
        try:
            handle = open(self.filename, 'rU')
        except IOError as e:
            sys.stderr.write("IOError: " + str(e) +
                             "\nError in SlnFile.parse_sln_file()\n")
            sys.stderr.flush()
            sys.exit()
        sln_list = []
        next_line = False
        for line in handle:
            line = line.strip('\n')
            if re.search(".txt", line):
                next_line = True
                continue
            if next_line:
                sln_list.append(line)
                next_line = False
        handle.close()
        return sln_list

    def sln_to_mol(self):
        '''
        Converts string-based sln

        :param sln_list: sln in string format in list
        :return: rdkit mol class in list
        '''
        mol_list = []
        for sln in self.sln_list:
            try:
                mol = MolFromSLN(sln)
            except ValueError:
                # ToDo Error tracking output at some point
                continue
            mol_list.append(mol)
        return filter(lambda x: x != None, mol_list)

    def mol_to_fp(self):
        if not self.mol_list:
            self.set_mol_list(self.sln_to_mol())
        fp_list = map(lambda x:
                      AllChem.GetMorganFingerprintAsBitVect(x, 2),
                      self.mol_list)
        #fp_list = map(lambda x: FingerprintMols.FingerprintMol(x),
        #              self.mol_list)
        return filter(lambda x: x != None, fp_list)

    def mol_to_plain_fp(self):
        if not self.mol_list:
            self.set_mol_list(self.sln_to_mol())
        fp_list = map(lambda x:
                      AllChem.GetMorganFingerprint(x, 2), self.mol_list)
        return filter(lambda x: x != None, fp_list)

    # Setters

    def set_mol_list(self, mol_list):
        self.mol_list = mol_list

    def set_fp_list(self, fp_list):
        self.fingerprint_list = fp_list

    # Getters

    def get_mol_list(self):
        if not self.mol_list:
            self.set_mol_list(self.sln_to_mol())
        return self.mol_list

    def get_fingerprint_list(self):
        if not self.fingerprint_list:
            self.set_fp_list(self.mol_to_fp())
        return self.fingerprint_list

    def get_plain_fingerprint_list(self):
        return self.mol_to_plain_fp()


class SmilesFile():
    def __init__(self, filename):
        self.filename = filename
        smile_list = []
        try:
            handle = open(self.filename, 'rU')
        except IOError as e:
            sys.stderr.write("IOError: " + str(e) +
                             "\nError in FileHandler.__init__()\n")
            sys.exit()
        for line in handle:
            line = line.strip('\n').split('\t')
            smile_str = line[3]
            smile_list.append(smile_str)
        handle.close()
        self.smile_list = smile_list
        self.mol_list = None
        self.fp_list = None

    def smile_to_mol(self):
        mol_list = []
        for smile in self.smile_list:
            try:
                mol = MolFromSmiles(smile)
            except:
                continue
            mol_list.append(mol)
        if len(mol_list) == 0:
            mol_list = None
        return mol_list


    def mol_to_fp(self):
        fp_list = []
        if not self.mol_list:
            self.set_mol_list(self.smile_to_mol())
        for mol in self.mol_list:
            try:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
                fp_list.append(fp)
            except:
                continue
        return fp_list


    # Setters

    def set_mol_list(self, mol_list):
        self.mol_list = mol_list


    def set_fp_list(self, fp_list):
        self.fp_list = fp_list


    # Getters

    def get_smile_list(self):
        return self.smile_list


    def get_mol_list(self):
        if not self.mol_list:
            self.set_mol_list(self.smile_to_mol())
        return self.mol_list


    def get_fingerprint_list(self):
        if not self.fp_list:
            self.set_fp_list(self.mol_to_fp())
        return self.fp_list
