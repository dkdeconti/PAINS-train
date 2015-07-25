__author__ = 'ddeconti'

import pickle
import FileHandler
import sys
import numpy
import random

from bokeh.plotting import figure, output_file, show, VBox, HBox

from rdkit import DataStructs
from rdkit.Chem import AllChem, SDMolSupplier

def randomly_pick_from_sdf(sdf_filename):
    sdf_struct = SDMolSupplier(sdf_filename)
    print len(sdf_struct)
    sdf_struct = random.sample(sdf_struct, 4000)
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


def fill_array(fps):
    np_fps = []
    for fp in fps:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    return np_fps


def fill_points(target, control, rf):
    x = []
    y = []
    for i in range(0, 101, 5):
        sys.stdout.write('.')
        sys.stdout.flush()
        p = len(target)
        n = len(control)
        fp = 0
        fn = 0
        tp = 0
        tn = 0
        for test in target:
            if rf.predict_proba(test)[0][1]*100 > i:
                tp += 1
            else:
                fn += 1
        for test in control:
            if rf.predict_proba(test)[0][1]*100 > i:
                fp += 1
            else:
                tn += 1
        specificity = tn/float(tn+fp)
        fpr = 1-specificity
        sensitivity = tp/float(tp+fn)
        x.append(fpr)
        y.append(sensitivity)
    sys.stdout.write('\n')
    return x, y


def main(sa):
    rf_filename = sa[0]
    pains_filename = sa[1]
    control_filename = sa[2]
    rf = pickle.load(open(rf_filename, 'rb'))
    pains_fps = fill_array(FileHandler.SlnFile(pains_filename)
                           .get_fingerprint_list())
    control_fps = fill_array(randomly_pick_from_sdf(control_filename))
    x, y = fill_points(random.sample(pains_fps, 40),
                       random.sample(control_fps, 400),
                       rf)
    output_file("rf_roc.html")
    p = figure(x_axis_label="False Positive Rate",
               y_axis_label="True Positive Rate")
    p.line(x, y)
    show(p)
    pass


if __name__ == "__main__":
    main(sys.argv[1:])
