__author__ = 'ddeconti'

import FileHandler
import numpy
import random
import pickle
import sys
from rdkit import DataStructs
from rdkit.Chem import AllChem, SDMolSupplier
from sklearn.cross_validation import train_test_split
from sklearn.neighbors import KNeighborsClassifier


def optimize_knn(target_train, target_test, control_train, control_test):
    out_list = ["num_trees", "specificity", "sensitivity", "fdr", "acc", "f1",
                "precision"]
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
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(test, arr)
        test = arr
        if knn.predict(test) == 1:
            tp += 1
        else:
            fn += 1
    for test in control:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(test, arr)
        test = arr
        if knn.predict(test) == 1:
            fp += 1
        else:
            tn += 1
    specificity = tn/float(tn+fp)
    sensitivity = tp/float(tp+fn)
    fdr = fp/float(tp+fp)
    acc = (tp+tn)/float(p+n)
    f1 = (2*tp)/float(2*tp+fp+fn)
    precision = (tp)/float(tp+fp)
    out_list = [specificity, sensitivity, fdr, acc, f1, precision]
    return out_list


def randomly_pick_from_sdf(sdf_filename, max_N=4000):
    sdf_struct = SDMolSupplier(sdf_filename)
    print len(sdf_struct)
    sdf_struct = random.sample(sdf_struct, max_N)
    try:
        mol_list = [m for m in sdf_struct]
    except:
        sys.stderr.write("Error parsing SDMolSupplier object\n" +
                         "Error in randomly_pick_from_sdf()\n")
        sys.exit()
    fp_list = []
    for m in mol_list:
        try:
            fp_list.append(AllChem.GetMorganFingerprintAsBitVect(m, 2))
        except:
            continue
    return filter(lambda x: x != None, fp_list)


def main(sa):
    sln_filename = sa[0]
    sdf_filename = sa[1]
    test_filename = sa[2]
    sln_fp = FileHandler.SlnFile(sln_filename).get_fingerprint_list()
    sdf_fp = randomly_pick_from_sdf(sdf_filename, max_N=4000)
    test_fp = FileHandler.SdfFile(test_filename).get_fingerprint_list()
    pain_train, pain_test = train_test_split(sln_fp,
                                             test_size=0.2,
                                             random_state=24)
    control_train, control_test = train_test_split(sdf_fp,
                                                   test_size=0.2,
                                                   random_state=24)
    #optimize_knn(pain_train, pain_test, control_train, control_test)
    knn = train_knn(pain_train, control_train, 5)
    stat_list = test_knn(sln_fp, test_fp, knn)
    out_list = ["num_trees", "specificity", "sensitivity", "fdr", "acc", "f1",
                "precision"]
    print '\t'.join(out_list)
    print '\t'.join(str(j) for j in stat_list)


if __name__ == "__main__":
    main(sys.argv[1:])