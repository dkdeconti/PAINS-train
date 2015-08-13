__author__ = 'ddeconti'

import FileHandler
import numpy
import pickle
import sys

from rdkit.Chem import DataStructs


def test_rf(target, control, rf):
    p = len(target)
    n = len(control)
    fp = 0
    fn = 0
    tp = 0
    tn = 0
    for test in target:
        if rf.predict(test) == 1:
            tp += 1
        else:
            fn += 1
    for test in control:
        if rf.predict(test) == 1:
            fp += 1
        else:
            tn += 1
    specificity = tn/float(tn+fp)
    sensitivity = tp/float(tp+fn)
    fdr = fp/float(tp+fp)
    acc = (tp+tn)/float(p+n)
    f1 = (2*tp)/float(2*tp+fp+fn)
    out_list = [specificity, sensitivity, fdr, acc, f1]
    return out_list


def fill_array(fps):
    np_fps = []
    for fp in fps:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    return np_fps


def main(sa):
    rf_filename = sa[0]
    pains_filename = sa[1]
    test_filename = sa[2]
    rf = pickle.load(open(rf_filename, "rb"))
    pains_fps = FileHandler.SlnFile(pains_filename).get_fingerprint_list()
    test_fps = FileHandler.SdfFile(test_filename).get_fingerprint_list()
    pains_fps = fill_array(pains_fps)
    test_fps = fill_array(test_fps)

    stat_list = test_rf(fill_array(pains_fps),
                        fill_array(test_fps),
                        rf)
    print stat_list


if __name__ == "__main__":
    main(sys.argv[1:])