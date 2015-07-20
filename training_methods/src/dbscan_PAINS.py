__author__ = 'ddeconti'

import FileHandler
import pickle
import numpy
import numpy as np
import sys

from sklearn.cluster import DBSCAN
from sklearn import metrics



def dbscan_array(pains):

    db = DBSCAN(eps=0.3, min_samples=10).fit(pains)
    core_samples_mask = numpy.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    print db.labels_

    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    import matplotlib.pyplot as plt
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    for k, col in zip(unique_labels, colors):
        if k == -1:
            col = 'k'
        class_member_mask = (labels == k)
        xy = pains[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=14)
        xy = pains[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                markeredgecolor='k', markersize=6)
    plt.show()
    pass


def main(sa):
    sln_filename = sa[0]
    try:
        sln_file = pickle.load(open("sln_file.p", "rb"))
    except:
        sln_file = FileHandler.SlnFile(sln_filename)
        pickle.dump(sln_file, open("sln_file.p", "wb"))
    try:
        sln_mol = pickle.load(open("sln_mol.p", "rb"))
    except:
        sln_mol = sln_file.get_mol_list()
        pickle.dump(sln_mol, open("sln_mol.p", "wb"))
    sln_fp = sln_file.get_fingerprint_list()
    dbscan_array(sln_fp)
    pass


if __name__ == "__main__":
    main(sys.argv[1:])