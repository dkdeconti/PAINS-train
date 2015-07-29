__author__ = 'ddeconti'


import FileHandler
import pandas
import random
import sys

from bokeh.palettes import Blues9
from bokeh.charts import HeatMap, output_file, show
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem, DataStructs, SDMolSupplier, Fingerprints


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
    return filter(lambda x: x != None, mol_list)


def mol_to_fp(mol_list):
    fp_list = []
    for mol in mol_list:
        try:
            fp = FingerprintMols.FingerprintMol(mol)
            fp_list.append(fp)
        except:
            continue
    fp_list = filter(lambda x: x != None, fp_list)
    return fp_list


def make_sim_matrix(fps):
    m = []
    for i in xrange(len(fps)):
        temp = []
        for j in xrange(len(fps)):
            temp.append(DataStructs.FingerprintSimilarity(fps[i], fps[j]))
        m.append(temp)
    return m

def plot_heatmap(m):
    import matplotlib.pyplot as plt
    import numpy as np
    m = [np.array(i) for i in m]
    m = np.array(m)
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(m)
    plt.savefig("heatmap.png")

def b_heatmap(df):
    output_file("heatmap.html")
    p = HeatMap(df, palette=Blues9)
    show(p)


def matrix_to_pd(m):
    d = {}
    idx = []
    for i in xrange(len(m)):
        idx.append(i)
        d[i] = m[i]
    df = pandas.DataFrame(d, index=idx)
    return df


def main(sa):
    sln_filename = sa[0]
    sdf_filename = sa[1]
    sln_mols = FileHandler.SlnFile(sln_filename).get_mol_list()
    sdf_mols = randomly_pick_from_sdf(sdf_filename, max_N=400)
    sln_fps = mol_to_fp(sln_mols)
    sdf_fps = mol_to_fp(sdf_mols)
    fps = sln_fps + sdf_fps
    sims = make_sim_matrix(fps)
    plot_heatmap(sims)
    #b_heatmap(matrix_to_pd(sims))


if __name__ == "__main__":
    main(sys.argv[1:])
