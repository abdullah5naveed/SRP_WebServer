import os
import csv
import numpy as np
from sklearn import decomposition
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from flask import Flask, render_template, request
from wtforms import Form, TextAreaField, validators
import extractFeatures
from keras.models import load_model
import joblib
import pickle
import time

project_root = os.path.dirname(os.path.realpath('__file__'))
template_path = os.path.join(project_root, 'templates')
static_path = os.path.join(project_root, 'static')
app = Flask(__name__, template_folder=template_path, static_folder=static_path)

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

class PredForm(Form):
    sequence = TextAreaField(u'\n\rEnter Sequence:', [validators.DataRequired()])

def SimpleParser(sequence):
    seq = sequence.split('\n')
    re = ''
    for x in seq:
        re = re + x[:len(x) - 2]
    return re

def SimpleFastaParser(fasta_sequence):
    seq = fasta_sequence.split('\n')
    seq = seq[1:]
    re = ''
    for x in seq:
        re = re + x[:len(x) - 2]
    return re

@app.route("/server", methods=['GET', 'POST'])
def index():
    form = PredForm(request.form)
    print(form.errors)
    if request.method == 'POST':
        input_seq = request.form['sequence']
        
        results = []
        
        seqs = input_seq.split('>')
        # loop  here
        model = joblib.load(('model.pkl'))
        model.set_params(n_jobs = 1)

        inputSize = 153
        dataset = np.genfromtxt("./FVs.csv", delimiter=",", dtype=float)
        X = np.array(dataset[:, 0:inputSize], dtype=np.float32)
        X = np.nan_to_num(X)
        std_scale = StandardScaler().fit(X)
        X = X.reshape(-1,153)
        X = std_scale.transform(X)
        X = np.array(X, dtype=np.float32)
        X = np.nan_to_num(X)
        pca = decomposition.PCA(n_components=2)
        pca.fit(X)
        X = pca.transform(X)
#            model = pickle.load(open('model.pkl','rb'))
        time.sleep(5)

        for ss in seqs[1:]:
            ss = '>'+ss
            sequence = SimpleFastaParser(ss)
            #sequence = SimpleParser(input_seq)
            featur = extractFeatures.get_features([sequence])
            np.random.seed(5)
            featur = np.array(featur, dtype=np.float32)
            featur = np.nan_to_num(featur)
            featur = std_scale.transform(featur)
            featur = np.nan_to_num(featur)
            featur = pca.transform(featur)
            r = model.predict(featur)
            sc = model.predict_proba(featur)
        
            if r[0] == 1.0:
                class1 = 'STRESS RESPONSE PROTEINS Sequences'
            if r[0] == 0.0:
                class1 = 'Non STRESS RESPONSE PROTEINS Sequences'

            result = [sequence, class1,np.max(sc)]
            results.append(result)
        return resultPage(results)

    return render_template('server.html', form=form, title="server")


def resultPage(result):
    return render_template('result.html', result=result, title="Results")


@app.route('/')
def intro():
    return render_template('intro.html', title="intro")



@app.route('/srp')
def srp():
    return render_template('srp.html', title="srp")


@app.route('/sample')
def sample():
    return render_template('sample.html', title="sample")

if __name__ == "__main__":
    app.run()
