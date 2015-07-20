__author__ = 'ddeconti'

import FileHandler
import pickle
import sys
from bokeh.palettes import Blues9
from bokeh.plotting import output_file, figure, show, VBox, HBox
from bokeh.charts import Histogram, HeatMap
from rdkit import DataStructs


def similarity_compare(fp):
    tanimoto_matrix = [[1] * len(fp)] * len(fp)
    for i in xrange(len(fp)):
        for j in xrange(len(fp)):
            if i == j:
                continue
            sim = DataStructs.FingerprintSimilarity(fp[i],
                                                    fp[j])
            tanimoto_matrix[i][j] = sim
    return tanimoto_matrix


def get_similarities_list(m):
    sim_list = []
    for i in xrange(len(m)):
        if i >= len(m) - 1:
            continue
        for j in xrange(i + 1, len(m)):
            sim_list.append(m[i][j])
    return sim_list


def plot_histogram(pains, control):
    '''
    distributions = OrderedDict(PAINs=pains, Control=control)
    df = pandas.DataFrame(distributions)
    distributions = df.to_dict()
    for k, v in distributions.items():
        distributions[k] = v.values()
    '''
    df = {"PAINs": pains, "Control": control}
    output_file("histograms.html")
    hist = Histogram(df, bins=20, legend=True)
    return hist


def plot_heatmap(all):
    p = HeatMap(all, palette=Blues9)
    return p


def main(sa):
    sln_filename = sa[0] # PAINs
    sdf_filename = sa[1] # Control set
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
        pains_fp = pickle.load(open("pains_fp.p", "rb"))
        control_fp = pickle.load(open("control_fp.p", "rb"))
    except:
        pains_fp = sln_file.get_fingerprint_list()
        control_fp = sdf_file.get_fingerprint_list()
        pickle.dump(pains_fp, open("pains_fp.p", "wb"))
        pickle.dump(control_fp, open("control_fp.p", "wb"))
    sys.stdout.write("Tanimoto similarity of PAINs.\n")
    sys.stdout.flush()
    try:
        pains_tanimoto = pickle.load(open("pains_tanimoto.p", "rb"))
    except:
        pains_tanimoto = similarity_compare(pains_fp)
        pickle.dump(pains_tanimoto, open("pains_tanimoto.p", "wb"))
    sys.stdout.write("Tanimoto similarity of Control.\n")
    sys.stdout.flush()
    try:
        control_tanimoto = pickle.load(open("control_tanimoto.p", "rb"))
    except:
        control_tanimoto = similarity_compare(control_fp)
        pickle.dump(control_tanimoto, open("control_tanimoto.p", "wb"))
    sys.stdout.write("Tanimoto similarity of both.\n")
    sys.stdout.flush()
    try:
        all_tanimoto = pickle.load(open("all_tanimoto.p", "rb"))
    except:
        all_tanimoto = similarity_compare(pains_fp + control_fp)
        pickle.dump(all_tanimoto, open("all_tanimoto.p", "wb"))
    sys.stdout.write("Plotting histograms.\n")
    sys.stdout.flush()
    hist = plot_histogram(get_similarities_list(pains_tanimoto),
                          get_similarities_list(control_tanimoto))
    sys.stdout.write("Plotting heatmap\n")
    sys.stdout.flush()
    heatmap = plot_heatmap(all_tanimoto)
    output_file("Pains_vs_Control_plots.html")
    VBox(hist, heatmap)
    show()


if __name__ == "__main__":
    main(sys.argv[1:])