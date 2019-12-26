import numpy as np
from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix
from Bucket import Bucket
from convert import convert
from scipy.signal import hilbert
from scipy.ndimage import filters
from plot_confusion import plot_confusion_matrix
import scipy.fftpack
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import kurtosis, skew, moment
from sklearn.externals import joblib
import matplotlib.pyplot as plt

#Reading data
d0 = np.load('../bench_train/60Hz/Multiplos/bc0_0.npy')
d50 = np.load('../bench_train/60Hz/Multiplos/bc50_0.npy')
d30 = np.load('../bench_train/60Hz/Multiplos/bc30_0.npy')
d60 = np.load('../bench_train/60Hz/Multiplos/bc60_0.npy')
d70 = np.load('../bench_train/60Hz/Multiplos/bc70_0.npy')

#Join data
dados = np.concatenate((d0, d60, d70, d50, d30))

#Creating desired responses
b = Bucket(5*500, 0, 60, False)
y = b.createLabelsSVM()

############################################################
#Feature Extractor
X = convert(dados)

#Spliting data
X_train, X_test, y_train, y_test = train_test_split(
     X, y, test_size=0.3, random_state=1)

#Training RF Classifier
clf = RandomForestClassifier(random_state=1, 
                             criterion='entropy',
                             n_estimators=100,
                             max_features='auto')
clf.fit(X_train, y_train)

################################### PLOT CONFUSION MATRIX

y_pred = clf.predict(X_test)

cnf_matrix = confusion_matrix(y_test, y_pred)
np.set_printoptions(precision=2)

# Plot non-normalized confusion matrix
plt.figure(figsize=(10,10))
plot_confusion_matrix(cnf_matrix, classes=['1','2','3','4','5','6',
                                           '7','8','9','10','11','12','13',
                                           '14','15','16','17','18','19', 
                                           '20','21','22'],
                      normalize=True,
                      title='Confusion Matrix')

plt.show()

#Saving classifier
joblib.dump(clf, 'RF_Com_Ruido1.pkl')
