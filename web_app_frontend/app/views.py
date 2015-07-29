__author__ = 'ddeconti'

import base64
import numpy
import sys
import pickle
from io import BytesIO
from flask import Flask, render_template, request, make_response, send_file
from rdkit import DataStructs
from rdkit.Chem import AllChem, MolFromSmiles, Draw
#from app import app


app = Flask(__name__)



@app.route('/')
@app.route('/index')
def index():
    return render_template("input.html")

@app.route('/about')
def pains_train_about():
    return render_template("about.html")

@app.route('/contact')
def pains_train_contact():
    return render_template("contact.html")

@app.route('/input')
def smile_input():
    return render_template("input.html")


@app.route('/output')
def smile_output():
    smile_str = request.args.get("smile_str")
    smile_str = str(smile_str)
    rf = pickle.load(open("static/rf_n300.p", "rb"))
    try:
        mol = MolFromSmiles(smile_str)
    except:
        return render_template("error.html",
                               error_str="Error in MolFromSmiles")
    sys.stdout.write(str(mol) + "\n")
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    except:
        return render_template("error.html",
                               error_str="Error in AllChem.GetMorgan" +
                                         "FingerprintAsBitVect()")
    sys.stdout.write(str(fp) + '\n')
    try:
        arr = numpy.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        fpm = arr
    except:
        return render_template("error.html",
                               error_str="Error in DataStructs." +
                               "ConvertToNumpyArray()")
    figfile = BytesIO()
    Draw.MolToFile(mol, figfile, imageType="png")
    figfile.seek(0)
    figdata_png = figfile.getvalue()
    figdata_png = base64.b64encode(figdata_png)
    p = rf.predict(fpm)
    prob = rf.predict_proba(fpm)[0][1]
    percent = int(prob*100)
    if p == 1:
        outstr = "promiscuous molecule"
    else:
        outstr = "non-promiscuous molecule"
    sys.stdout.write(smile_str + "\n")
    return render_template("output.html", smile_str=smile_str,
                           binary_str=outstr, img_obj=figdata_png,
                           predict_prob=percent)



@app.route('/error')
def error_output():
    return render_template("error.html")


if __name__ == "__main__":
    app.run(host='0.0.0.0', port=5000, debug=True)
