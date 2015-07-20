__author__ = 'ddeconti'

import FileHandler
import pickle
import numpy as np
import sys

from sklearn import cluster
from rdkit import DataStructs
from matplotlib import pyplot

def plot_kmeans(pains):
    np_fps = []
    for fp in pains:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    k = 2
    kmeans = cluster.KMeans(n_clusters=k)
    kmeans.fit(np_fps)
    labels = kmeans.labels_
    centroids = kmeans.cluster_centers_
    for i in range(k):
        #print i, np.where(labels==i)
        #print np.where(labels==i)[0]
        ds = np.array([np_fps[j] for j in np.where(labels==i)[0]])
        #ds = np_fps[np.where(labels==i)[0]]
        print ds
        pyplot.plot(ds[:, 0], ds[:, 1], 'o')
        lines = pyplot.plot(centroids[i,0], centroids[i,1], 'kx')
        pyplot.setp(lines,ms=15.0)
        pyplot.setp(lines,mew=2.0)
    pyplot.show()

def main(sa):
    sln_filename = sa[0]
    try:
        sln_file = pickle.load(open("sln_file.p", "rb"))
    except:
        sln_file = FileHandler.SlnFile(sln_filename)
        pickle.dump(sln_file, open("sln_file.p", "wb"))
    sln_fp = sln_file.get_fingerprint_list()
    plot_kmeans(sln_fp)


if __name__ == "__main__":
    main(sys.argv[1:])