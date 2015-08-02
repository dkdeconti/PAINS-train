__author__ = 'ddeconti'


import FileHandler
import numpy
import sys

from bokeh.plotting import figure, output_file, show, VBox, HBox
from rdkit import DataStructs
from sklearn.cluster import KMeans
from sklearn.decomposition.pca import PCA
from sklearn.metrics import silhouette_score


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


def train_kmeans(pc_array, n_clusters):
    k_means = KMeans(n_clusters=n_clusters)
    k_means.fit(pc_array[:,0:2])
    return k_means


def plot_k_means(k_means, pc_array):
    plot_dict = {}
    for i in xrange(len(k_means.labels_)):
        kl = k_means.labels_[i]
        if kl in plot_dict:
            plot_dict[kl]["PC1"].append(pc_array[i][0])
            plot_dict[kl]["PC2"].append(pc_array[i][1])
        else:
            plot_dict[kl] = {}
            plot_dict[kl]["PC1"] = [pc_array[i][0]]
            plot_dict[kl]["PC2"] = [pc_array[i][1]]
    output_file("pca_kmeans_cluster.html")
    p = figure(x_axis_label = "PC1",
               y_axis_label = "PC2")
    color_vector = ["blue", "red", "green", "purple", "orange", "cyan",
                    "magenta", "yellow", "black"]
    for k in sorted(plot_dict.keys()):
        p.scatter(plot_dict[k]["PC1"], plot_dict[k]["PC2"],
                  color = color_vector[k])
    return p


def recurr_plots(reduced_pains_fps, current_iter, max_iter=9):
    kmeans_predict = train_kmeans(reduced_pains_fps, current_iter)
    p = plot_k_means(kmeans_predict, reduced_pains_fps)
    if current_iter == max_iter:
        return p
    else:
        current_iter += 1
        return VBox(p, recurr_plots(reduced_pains_fps, current_iter))


def kmeans_sil_analysis(reduced_pains_fps):
    s = []
    for n_clusters in range(2, 21):
        kmeans = train_kmeans(reduced_pains_fps, n_clusters)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        s.append(silhouette_score(reduced_pains_fps[:,:],
                                  labels,
                                  metric="euclidean"))
    print s
    p = figure(title="Silhouette Scoring of PCA K-means clusters",
               x_axis_label="Number of k clusters",
               y_axis_label="Silhouette score")
    p.line(range(2,21), s)
    return p


def main(sa):
    pains_filename = sa[0]
    sdf_filename = sa[1]
    pains_fps = FileHandler.SlnFile(pains_filename).get_fingerprint_list()
    sdf_fps = FileHandler.SdfFile(sdf_filename).get_fingerprint_list()

    reduced_pains_fps = train_pca(pains_fps+sdf_fps)
    #kmeans_predict = train_kmeans(reduced_pains_fps)
    #plot_k_means(kmeans_predict, reduced_pains_fps)

    plots = recurr_plots(reduced_pains_fps, 1)
    sil_plot = kmeans_sil_analysis(reduced_pains_fps)
    show(VBox(plots, sil_plot))


if __name__ == "__main__":
    main(sys.argv[1:])
