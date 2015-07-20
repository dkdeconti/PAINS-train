__author__ = 'ddeconti'

import FileHandler
import numpy
import pickle
import sys
from bokeh.plotting import figure, output_file, show, VBox, HBox
from rdkit import DataStructs
from sklearn.decomposition.pca import PCA
from sklearn.cross_validation import train_test_split


def pca(target, control, title):
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
               color="blue")
    p1.scatter(np_fps_r[len(target):, 0], np_fps_r[len(target):, 1],
               color="red")
    p2 = figure(x_axis_label="PC2",
                y_axis_label="PC3",
                title=title)
    p2.scatter(np_fps_r[:len(target), 1], np_fps_r[:len(target), 2],
               color="blue")
    p2.scatter(np_fps_r[len(target):, 1], np_fps_r[len(target):, 2],
               color="red")
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


def main(sa):
    sln_filename = sa[0]
    sdf_filename = sa[1]
    smiles_filename = sa[2]
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
        smile_file = pickle.load(open("smile_file.p", 'rb'))
    except:
        smile_file = FileHandler.SmilesFile(smiles_filename)
        pickle.dump(smile_file, open("smile_file.p", "wb"))
    try:
        sln_mol = pickle.load(open("sln_mol.p", "rb"))
        sdf_mol = pickle.load(open("sdf_mol.p", "rb"))
        smile_mol = pickle.load(open("smile_mol.p", "rb"))
    except:
        sln_mol = sln_file.get_mol_list()
        sdf_mol = sdf_file.get_mol_list()
        smile_mol = sdf_file.get_mol_list()
        pickle.dump(sln_mol, open("sln_mol.p", "wb"))
        pickle.dump(sdf_mol, open("sdf_mol.p", "wb"))
        pickle.dump(smile_mol, open("smile_mol.p", "wb"))
    sln_fp = sln_file.get_fingerprint_list()
    sdf_fp = sdf_file.get_fingerprint_list()
    smile_fp = smile_file.get_fingerprint_list()

    '''
    pain_train, pain_test = train_test_split(sln_fp,
                                             test_size=0.2,
                                             random_state=24)
    control_train, control_test = train_test_split(sdf_fp,
                                                   test_size=0.2,
                                                   random_state=24)
    '''
    print "PCA for PAINS vs. Control"
    pvc = pca(sln_fp, sdf_fp, "PAINS vs. Screen Library")
    print "PCA for PAINS vs. All"
    pva = pca(sln_fp, smile_fp, "PAINS vs. Diverse Chemical Set")
    print "PCA within PAINS"
    pvp = pca_no_labels(sln_fp)
    output_file("pca_plots.html")
    p = VBox(pvc, pva, pvp)
    show(p)

if __name__ == "__main__":
    main(sys.argv[1:])