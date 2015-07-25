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
    sln_fp = FileHandler.SlnFile(sln_filename).get_fingerprint_list()
    sdf_fp = FileHandler.SdfFile(sdf_filename).get_fingerprint_list()
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