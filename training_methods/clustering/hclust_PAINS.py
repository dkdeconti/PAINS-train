__author__ = 'ddeconti'

import FileHandler
import matplotlib.pyplot as plt
import numpy
import sys
import scipy.cluster.hierarchy as hcluster

from bokeh.plotting import figure, output_file, show, VBox, HBox
from rdkit import DataStructs
from sklearn.decomposition.pca import PCA


def pca_plot(fp_list, clusters):
    np_fps = []
    for fp in fp_list:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    pca = PCA(n_components=3)
    pca.fit(np_fps)
    np_fps_r = pca.transform(np_fps)
    p1 = figure(x_axis_label="PC1",
                y_axis_label="PC2",
                title="PCA clustering of PAINS")
    p2 = figure(x_axis_label="PC2",
                y_axis_label="PC3",
                title="PCA clustering of PAINS")
    color_vector = ["blue", "red", "green", "orange", "pink", "cyan", "magenta",
                    "brown", "purple"]
    print len(set(clusters))
    for clust_num in set(clusters):
        print clust_num
        local_cluster = []
        for i in xrange(len(clusters)):
            if clusters[i] == clust_num:
                local_cluster.append(np_fps_r[i])
        print len(local_cluster)
        p1.scatter(np_fps_r[:,0], np_fps_r[:,1],
                   color=color_vector[clust_num])
        p2.scatter(np_fps_r[:,1], np_fps_r[:,2],
                   color=color_vector[clust_num])
    return HBox(p1, p2)


def clust(fp_list):
    np_fps = []
    for fp in fp_list:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    thresh = 6.5
    clusters = hcluster.fclusterdata(np_fps, thresh, criterion="distance")
    return clusters


def main(sa):
    pains_filename = sa[0]
    fp_list = FileHandler.SlnFile(pains_filename).get_fingerprint_list()
    clusters = clust(fp_list)
    p = pca_plot(fp_list, clusters)
    output_file("PCA_w_hclust.html")
    show(p)


if __name__ == "__main__":
    main(sys.argv[1:])