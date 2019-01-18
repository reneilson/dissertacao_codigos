import numpy as np
from sklearn.cross_validation import train_test_split
from sklearn.metrics import confusion_matrix
from Bucket import Bucket
from scipy.signal import hilbert
from scipy.ndimage import filters
from plot_confusion import plot_confusion_matrix
import scipy.fftpack
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import kurtosis, skew, moment
from sklearn.externals import joblib

#Lendo dados do SED-Bench
d30 = np.load('/bench_train/60Hz/Multiplos/bc30_0.npy')
d20 = np.load('/bench_train/60Hz/Multiplos/bc20_0.npy')

#Concatenando dados
dados = np.concatenate((d20, d30))

#Criando vetor de rotulos com codigo fonte do SED-Bench
b = Bucket(1000, 30,60, False)
y = b.createLabelsSVM()

#Deletando variaveis nao utilizaveis
del b, d20, d30

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

#Filtrando FFT dos sinais 
df2 = filters.gaussian_filter1d(X_,10)

#Calculando variaveis
eng_data = (dados**2).sum(axis=1)
max_df = np.max(abs(df), axis=1)
var_X2 = np.var(X_2[:,2:10],axis=1)
max_mean_X2 = np.max(X_[:,3:], axis=1)-np.mean(X_[:,3:], axis=1)
sum_X = np.sum((100*X_[:,3:])**4, axis=1)
X2_sumX2 = (X_2[:,2])*np.sum(X_2[:,:10], axis=1)
max_X = np.max(X_[:,3:], axis=1)
max_filter_X = np.max(df2[:,50:],axis=1)
sum_dados_df = np.sum((dados-df)**2, axis=1)
kurt_filter_X = kurtosis(df2, axis=1)
skew_filter_X = skew(df2, axis=1)
moment4_dados = moment(dados, 4, axis=1)
moment2_dados = moment(dados, 2, axis=1)
skew_std_sum_X = skew(X_[:,4:], axis=1)+1000*np.std(X_[:,3:],axis=1)+10*np.sum(X_[:,4:], axis=1)
max_filter_X_sum_X = np.max(df2[:,50:],axis=1)*np.sum((100*X_[:,3:])**4, axis=1)
first_dado_df = (dados-df)[:,0]
first_X_X_2 = X_[:,0]+X_2[:,0]
max_mean_X3 = np.max(X_3[:,4:], axis=1)-np.mean(X_3, axis=1)
kurt_hilbert = kurtosis(amplitude_envelope, axis=1)
moment_2_4_hilbert = moment(amplitude_envelope,  2, axis=1)*moment(amplitude_envelope,  4, axis=1)
moment_2_X2 = moment(X_2[:,3:], 2, axis=1)
max_min_mean_hilbert = (np.max(amplitude_envelope, axis=1) +np.min(amplitude_envelope, axis=1))*np.mean(amplitude_envelope, axis=1)
sum_filter_dados = np.sum(df, axis=1)
std_filter_X = np.std(df2, axis=1)

###########################################################
#Montando Vetor de Caracteristicas
X = np.array([max_df,
              eng_data,
              var_X2,
              max_mean_X2,
              sum_X,
              X2_sumX2,
              max_X,
              max_filter_X,
              sum_dados_df,
              kurt_filter_X,
              skew_filter_X,
              moment4_dados,
              moment2_dados,
              skew_std_sum_X,
              max_filter_X_sum_X,
              first_dado_df,
              first_X_X_2,
              max_mean_X3,
              kurt_hilbert,
              moment_2_4_hilbert,
              moment_2_X2,
              max_min_mean_hilbert,
              sum_filter_dados,
              std_filter_X,
              np.sum(X_[:,1:10], axis=1)*kurtosis(amplitude_envelope, axis=1),
              np.std(X_[:,1:10], axis=1)*np.var(dados, axis=1),
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
              dados[:,0],
              df[:,0],
              np.sum(df2, axis=1)*skew(X_, axis=1),
              np.std(dados-df, axis=1),
              (kurtosis(X_, axis=1) - kurtosis(df2, axis=1))*kurtosis(X_-df2, axis=1),
              ]).T

#Separando entre conjunto de treino e teste
X_train, X_test, y_train, y_test = train_test_split(
     X, y, test_size=0.3, random_state=1)

#Treinando Classificador ja com melhores hiperparametros
clf = RandomForestClassifier(random_state=1, criterion='entropy',
                             max_depth=15, n_estimators=50)
clf.fit(X_train, y_train)

#Salvando Classificador Treinado
joblib.dump(clf, 'RF_Com_Ruido.pkl')