__author__ = 'ddeconti'

import FileHandler
import numpy
import random
import sys
from bokeh.plotting import figure, output_file, show, VBox, HBox
from rdkit import DataStructs
from rdkit.Chem import AllChem, SDMolSupplier
from sklearn.decomposition.pca import PCA
from sklearn.cross_validation import train_test_split


def pca(target, control, title, name_one, name_two):
    np_fps = []
    for fp in target + control:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    ys_fit = [1] * len(target) + [0] * len(control)
    names = ["PAINS", "Control"]
    pca = PCA(n_components=3)
    pca.fit(np_fps)
    np_fps_r = pca.transform(np_fps)
    p1 = figure(x_axis_label="PC1",
                y_axis_label="PC2",
                title=title)
    p1.scatter(np_fps_r[:len(target), 0], np_fps_r[:len(target), 1],
               color="blue", legend=name_one)
    p1.scatter(np_fps_r[len(target):, 0], np_fps_r[len(target):, 1],
               color="red", legend=name_two)
    p2 = figure(x_axis_label="PC2",
                y_axis_label="PC3",
                title=title)
    p2.scatter(np_fps_r[:len(target), 1], np_fps_r[:len(target), 2],
               color="blue", legend=name_one)
    p2.scatter(np_fps_r[len(target):, 1], np_fps_r[len(target):, 2],
               color="red", legend=name_two)
    return HBox(p1, p2)

def pca_no_labels(target):
    np_fps = []
    for fp in target:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    pca = PCA(n_components=3)
    pca.fit(np_fps)
    np_fps_r = pca.transform(np_fps)
    p3 = figure(x_axis_label="PC1",
                y_axis_label="PC2",
                title="PCA clustering of PAINS")
    p3.scatter(np_fps_r[:, 0], np_fps_r[:, 1], color="blue")
    p4 = figure(x_axis_label="PC2",
                y_axis_label="PC3",
                title="PCA clustering of PAINS")
    p4.scatter(np_fps_r[:, 1], np_fps_r[:, 2], color="blue")
    return HBox(p3, p4)


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


def main(sa):
    sln_filename = sa[0]
    sdf_filename = sa[1]
    setsix_filename = sa[2]
    smiles_filename = sa[3]
    sln_fp = FileHandler.SlnFile(sln_filename).get_fingerprint_list()
    sdf_fp = randomly_pick_from_sdf(sdf_filename)
    setsix_fp = FileHandler.SdfFile(setsix_filename).get_fingerprint_list()
    smile_fp = FileHandler.SmilesFile(smiles_filename).get_fingerprint_list()
    print "PCA for PAINS vs. Chembl"
    pvc = pca(sln_fp, sdf_fp, "PAINS vs. ChEMBL",
              "PAINS", "ChEMBL")
    print "PCA for PAINS vs. Set six"
    pvb = pca(sln_fp, setsix_fp, "PAINS vs. ChEMBL set 5",
              "PAINS", "ChEMBL.5")
    print "PCA for PAINS vs. STitch"
    pva = pca(sln_fp, smile_fp, "PAINS vs. Stitch",
              "PAINS", "Stitch")
    print "PCA within PAINS"
    pvp = pca_no_labels(sln_fp)
    output_file("pca_plots.html")
    p = VBox(pvc, pvb, pva, pvp)
    show(p)

if __name__ == "__main__":
    main(sys.argv[1:])