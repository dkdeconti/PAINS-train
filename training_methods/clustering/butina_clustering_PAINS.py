__author__ = 'ddeconti'

import FileHandler
import sys

from rdkit import DataStructs
from rdkit.ML.Cluster import Butina


def cluster_fps(fps, cutoff=0.2):
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1-x for x in sims])
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


def main(sa):
    pains_filename = sa[0]
    fps = FileHandler.SlnFile(pains_filename).get_plain_fingerprint_list()
    cs = cluster_fps(fps, cutoff=0.4)
    print cs
    print len(cs), len(fps)
    pass


if __name__ == "__main__":
    main(sys.argv[1:])