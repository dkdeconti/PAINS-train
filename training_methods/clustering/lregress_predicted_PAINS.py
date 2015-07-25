__author__ = 'ddeconti'

import FileHandler
import math
import numpy
import pickle
import sys

from bokeh.plotting import figure, output_file, show, HBox
from rdkit.Chem import AllChem, MolFromSmiles


def plot_scatter(x, y, title, x_label, y_label, color="red"):
    plot = figure(x_axis_label=x_label,
                  y_axis_label=y_label)
    plot.scatter(x, y, fill_color="red")
    return plot


def build_prediction_result_array(rf, drug_list):
    '''
    Attemps loading pickle (may need to turn off function

    :param rf:
    :param filename:
    :return:
    '''
    prediction_probs = []
    for drug in drug_list:
        try:
            mol = MolFromSmiles(drug.get_SMILES())
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
        except:
            continue
        p = rf.predict_proba(fp)[0][1]
        drug.set_rf_prediction_prob(p)
    return drug_list


def polyfit(x, y, degree):
    new_x = []
    new_y = []
    for i in xrange(len(x)):
        if x[i] != None and y[i] != None:
            new_x.append(x[i])
            new_y.append(y[i])
    x = new_x
    y = new_y
    results = {}
    coeffs = numpy.polyfit(x, y, degree)
    results['polynomial'] = coeffs.tolist()
    correlation = numpy.corrcoef(x,y)[0,1]
    results['correlation'] = correlation
    results['determination'] = correlation**2
    return results


def main(sa):
    rf_pickle_name = sa[0]
    chem_dict_filename = sa[1]
    try:
        rf = pickle.load(open(rf_pickle_name, "rb"))
    except IOError as e:
        sys.stderr.write("IOError: " + str(e) +
                         "\nError in main()\n")
        sys.exit()
    drug_list = FileHandler.WikiScrapedDB(chem_dict_filename).get_drug_list()
    drug_list = build_prediction_result_array(rf, drug_list)
    p1 = plot_scatter([drug.get_rf_prediction_prob() for drug in drug_list],
                      [drug.get_side_effect_count() for drug in drug_list],
                      "Chemical promiscuity vs. # side effects",
                      "Predicted probability of PAINS classification",
                      "Number of side effects")
    p2 = plot_scatter([drug.get_rf_prediction_prob() for drug in drug_list],
                      [drug.get_gene_effect_count() for drug in drug_list],
                      "Chemical promiscuity vs. # gene expression changes",
                      "Predicted probability of PAINS classification",
                      "Number of identified gene expression effects")
    results = polyfit([drug.get_rf_prediction_prob() for drug in drug_list],
                      [drug.get_side_effect_count() for drug in drug_list],
                      1)
    print results
    output_file("lregress.html")
    show(HBox(p1, p2))

if __name__ == "__main__":
    main(sys.argv[1:])