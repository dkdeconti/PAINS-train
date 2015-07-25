__author__ = 'ddeconti'


import FileHandler
import numpy
import sys

from bokeh.plotting import figure, output_file, show, VBox, HBox
from rdkit import DataStructs
from sklearn.cluster import DBSCAN
from sklearn.decomposition.pca import PCA


def train_pca(pains_fps, num_components=3):
    '''
    Dimensional reduction of fps bit vectors to principal components
    :param pains_fps:
    :return: pca reduced fingerprints bit vectors
    '''
    np_fps = []
    for fp in pains_fps:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    pca = PCA(n_components=num_components)
    pca.fit(np_fps)
    fps_reduced = pca.transform(np_fps)
    return fps_reduced


def train_dbscan(pains_fps):
    db = DBSCAN(eps=1, min_samples=10).fit(pains_fps)
    print db.labels_
    pass


def main(sa):
    pains_filename = sa[0]
    pains_fps = FileHandler.SlnFile(pains_filename).get_fingerprint_list()
    reduced_pains_fps = train_pca(pains_fps, num_components=2)
    train_dbscan(reduced_pains_fps)
    p = figure(x_axis_label="PC1",
               y_axis_label="PC2")


if __name__ == "__main__":
    main(sys.argv[1:])
