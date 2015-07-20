__author__ = 'ddeconti'

import FileHandler
import numpy
import pickle
import sys
from rdkit import Chem, DataStructs
from sklearn import svm, datasets
from sklearn.cross_validation import train_test_split



def svc_training(target, control):
    np_fps = []
    for fp in target + control:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    lin_svc = svm.LinearSVC(C=1.0)
    ys_fit = [1] * len(target) + [0] * len(control)
    lin_svc.fit(np_fps, ys_fit)
    return lin_svc


def svc_test(target, control, lin_svc):
    ys_fit = [1] * len(target) + [0] * len(control)
    print lin_svc.score(target+control, ys_fit)


def main(sa):
    sln_filename = sa[0]
    sdf_filename = sa[1]
    try:
        sln_file = pickle.load(open("sln_file.p", "rb"))
    except:
        sln_file = FileHandler.SlnFile(sln_filename)
        pickle.dump(sln_file, open("sln_file.p", "wb"))
    try:
        sdf_file = pickle.load(open("sdf_file.p", "rb"))
    except:
        sdf_file = FileHandler.SdfFile(sdf_filename)
        pickle.dump(sdf_file, open("sdf_file.p", "wb"))
    try:
        sln_mol = pickle.load(open("sln_mol.p", "rb"))
        sdf_mol = pickle.load(open("sdf_mol.p", "rb"))
    except:
        sln_mol = sln_file.get_mol_list()
        sdf_mol = sdf_file.get_mol_list()
        pickle.dump(sln_mol, open("sln_mol.p", "wb"))
        pickle.dump(sdf_mol, open("sdf_mol.p", "wb"))
    sln_fp = sln_file.get_fingerprint_list()
    sdf_fp = sdf_file.get_fingerprint_list()
    pain_train, pain_test = train_test_split(sln_fp,
                                             test_size=0.2,
                                             random_state=24)
    control_train, control_test = train_test_split(sdf_fp,
                                                   test_size=0.2,
                                                   random_state=24)
    lin_svc = svc_training(pain_train, control_train)
    svc_test(pain_test, control_test, lin_svc)


if __name__ == "__main__":
    main(sys.argv[1:])