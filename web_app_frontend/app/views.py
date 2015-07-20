__author__ = 'ddeconti'

import base64
import pickle
from io import BytesIO
from flask import render_template, request, make_response, send_file
from rdkit.Chem import AllChem, MolFromSmiles, Draw
from rdkit.Chem.rdSLNParse import MolFromSLN
from app import app

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

@app.route('/')
@app.route('/index')
def index():
    return render_template("index.html")


@app.route('/input')
def smile_input():
    return render_template("input.html")


@app.route('/output')
def smile_output():
    smile_str = request.args.get("smile_str")
    rf = pickle.load(open("/home/ddeconti/PycharmProjects/PAINS-train/" +
                          "web_app_frontend/app/static/rf_n300.p", "rb"))
    try:
        mol = MolFromSmiles(smile_str)
    except:
        return render_template("error.html",
                               error_str="Error in MolFromSmiles")
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    except:
        return render_template("error.html",
                               error_str="Error in AllChem.GetMorgan" +
                                         "FingerprintAsBitVect()")
    figfile = BytesIO()
    Draw.MolToFile(mol, figfile, imageType="png")
    figfile.seek(0)
    figdata_png = figfile.getvalue()
    figdata_png = base64.b64encode(figdata_png)
    p = rf.predict(fp)
    if p == 1:
        outstr = "promiscuous molecule"
    else:
        outstr = "non-promiscuous molecule"
    return render_template("output.html", smile_str=smile_str,
                           binary_str=outstr, img_obj=figdata_png)


@app.route('/error')
def error_output():
    return render_template("error.html")



