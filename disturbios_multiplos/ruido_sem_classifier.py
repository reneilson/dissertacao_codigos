import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix
from Bucket import Bucket
from scipy.signal import hilbert
from scipy.ndimage import filters
import scipy.fftpack
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import kurtosis, skew, moment
from sklearn.externals import joblib
from Bucket import Bucket

#Lendo dados do SED-Bench
d0 = np.load('/bench_train/60Hz/Multiplos/bc0_0.npy')
d50 = np.load('/bench_train/60Hz/Multiplos/bc50_0.npy')

#Concatenando dados
dados = np.concatenate((d0, d50))

#Criando vetor de rotulos com codigo fonte do SED-Bench
b = Bucket(1000, 0, 60, False)
y = b.createLabelsSVM()

#Deletando variaveis nao utilizaveis
del b, d0, d50

#Filtrando dados de entrada
df=filters.gaussian_filter1d(dados,10)

#Calculando FFT dos sinais
Y = scipy.fftpack.fft(dados)
L = dados.shape[1]
P2 = np.abs(2*Y/L)
X_ = P2[:,1:np.round(L/2)+1]
del Y, L, P2

#Calculando FFT dos sinais filtrados
Y = np.fft.fft(df)
L = dados.shape[1]
P2 = np.abs(2*Y/L)
X_2 = P2[:,1:np.round(L/2)+1]
del Y, L, P2

#Calculando FFT dos sinais retirando o filtro
Y = np.fft.fft(dados-df)
L = (dados-df).shape[1]
P2 = np.abs(Y/L)
X_3 = P2[:,1:np.round(L/2)+1]
del Y, L, P2

#Calculando envelopamento dos sinais
analytic_signal = hilbert(dados)
amplitude_envelope = np.abs(analytic_signal)

#Calculando variaveis
eng_data = (dados**2).sum(axis=1)
max_df = np.max(abs(df), axis=1)
var_X2 = np.var(X_2[:,2:10],axis=1)
sum_dados_df = np.sum((dados-df)**2, axis=1)
moment4_dados = moment(dados, 4, axis=1)
moment2_dados = moment(dados, 2, axis=1)
first_dado_df = (dados-df)[:,0]
max_mean_X3 = np.max(X_3[:,4:], axis=1)-np.mean(X_3, axis=1)
kurt_hilbert = kurtosis(amplitude_envelope, axis=1)
sum_filter_dados = np.sum(df, axis=1)
std_filter_X = np.std(df, axis=1)

###########################################################
#Montando Vetor de Caracteristicas
X = np.array([max_df,
              eng_data,
              var_X2,
              sum_dados_df,
              moment4_dados,
              moment2_dados,
              first_dado_df,
              max_mean_X3,
              kurt_hilbert,
              sum_filter_dados,
              std_filter_X,
              np.sum(X_[:,1:10], axis=1)*kurtosis(amplitude_envelope, axis=1),
              np.sum(X_[:,10:], axis=1),
              np.std(X_[:,10:], axis=1),
              np.max(X_[:,10:], axis=1)-np.mean(X_[:,10:], axis=1),
              np.max(X_[:,50:100], axis=1)-np.mean(X_[:,50:100], axis=1),
              np.max(X_[:,150:200], axis=1)-np.mean(X_[:,150:200], axis=1),
              np.max(X_[:,200:], axis=1)-np.mean(X_[:,200:], axis=1),
              np.sum(X_[:,:50], axis=1)-
              np.sum(X_[:,50:100], axis=1)*np.sum(X_[:,150:200], axis=1),
              np.sum(X_[:,200:], axis=1),
              np.std(X_[:,:50], axis=1),
              np.std(X_[:,150:200], axis=1),
              X_[:,2],
              (dados-df)[:,0],
              dados[:,0],
              df[:,0],
              ]).T

#Separando entre conjunto de treino e teste
X_train, X_test, y_train, y_test = train_test_split(
         X, y, test_size=0.3, random_state=1)

#Treinando Classificador ja com melhores hiperparametros
clf = RandomForestClassifier(random_state=1, criterion='entropy',
                             max_depth=14, n_estimators=100)

clf.fit(X_train, y_train)

#Salvando Classificador Treinado
joblib.dump(clf, 'RF_Sem_Ruido.pkl')