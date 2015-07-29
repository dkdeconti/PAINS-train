__author__ = 'ddeconti'

import FileHandler
import random
import sys

from bokeh.charts import output_file, Histogram, show
from bokeh.models import Range1d
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


def get_same_similarity(fps):
    same_sims = []
    for i in xrange(len(fps)):
        for j in xrange(i+1, len(fps)):
            sim = DataStructs.FingerprintSimilarity(fps[i], fps[j])
            same_sims.append(sim)
    return same_sims


def get_diff_sims(pains, rest):
    diff_sims = []
    for i in xrange(len(pains)):
        for j in xrange(len(rest)):
            sim = DataStructs.FingerprintSimilarity(pains[i], rest[j])
            diff_sims.append(sim)
    return diff_sims


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


def plot_histogram(same, diff):
    df = {"PAINS vs. PAINS": same, "PAINS vs. ChEMBL": diff}
    output_file("histogram_pains_self_diff.html")
    hist = Histogram(df, bins=20, density=True, legend=True)
    hist.x_range = Range1d(0, 1)
    #hist.legend.orientation = "top_right"
    show(hist)


def main(sa):
    sln_filename = sa[0]
    sdf_filename = sa[1]
    sln_mols = FileHandler.SlnFile(sln_filename).get_mol_list()
    sdf_mols = randomly_pick_from_sdf(sdf_filename, max_N=400)
    print len(sln_mols), len(sdf_mols)
    sln_fps = mol_to_fp(sln_mols)
    sdf_fps = mol_to_fp(sdf_mols)
    same_sims = get_same_similarity(sln_fps)
    diff_sims = get_diff_sims(sln_fps, sdf_fps)
    plot_histogram(same_sims, diff_sims)


if __name__ == "__main__":
    main(sys.argv[1:])