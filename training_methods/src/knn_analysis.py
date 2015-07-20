__author__ = 'ddeconti'

import FileHandler
import numpy
import pickle
import sys
from rdkit import DataStructs
from sklearn.cross_validation import train_test_split
from sklearn.neighbors import KNeighborsClassifier


def optimize_knn(target_train, target_test, control_train, control_test):
    out_list = ["num_trees", "specificity", "sensitivity", "fdr", "acc", "f1"]
    print '\t'.join(out_list)
    for i in [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500,
              1000, 2000, 3000, 4000, 5000]:
        knn = train_knn(target_train, control_train, i)
        stat_list = test_knn(target_test, control_test, knn)
        print '\t'.join(str(j) for j in [i] + stat_list)


def train_knn(target, control, n_neigh):
    knn = KNeighborsClassifier()
    np_fps = []
    for fp in target + control:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    ys_fit = [1] * len(target) + [0] * len(control)
    knn.fit(np_fps, ys_fit)
    return knn


def test_knn(target, control, knn):
    p = len(target)
    n = len(control)
    fp = 0
    fn = 0
    tp = 0
    tn = 0
    ys_fit = [1] * len(target) + [0] * len(control)
    for test in target:
        if knn.predict(test) == 1:
            tp += 1
        else:
            fn += 1
    for test in control:
        if knn.predict(test) == 1:
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
    optimize_knn(pain_train, pain_test, control_train, control_test)
    #knn = train_knn(pain_train, control_train)
    #test_knn(pain_test, control_test, knn)


if __name__ == "__main__":
    main(sys.argv[1:])