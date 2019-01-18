import numpy as np
from sklearn.cross_validation import train_test_split
from scipy.ndimage import filters
import scipy.fftpack
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib

#Lendo dados do SED-Bench
d0 = np.load('/bench_train/60Hz/Multiplos/bc0_0.npy')
d50 = np.load('/bench_train/60Hz/Multiplos/bc50_0.npy')
d30 = np.load('/bench_train/60Hz/Multiplos/bc30_0.npy')
d20 = np.load('/bench_train/60Hz/Multiplos/bc20_0.npy')

#Concatenando dados
dados = np.concatenate((d0, d50, d30, d20))

#Criando vetor de rotulos 1-sem ruido 0-com ruido
y = np.concatenate((np.ones(500*2), np.zeros(500*2)))

#Deletando variaveis nao utilizaveis
del d0, d20, d30, d50

#Filtrando dados de entrada
df=filters.gaussian_filter1d(dados,10)

#Calculando FFT dos sinais
Y = scipy.fftpack.fft(dados)
L = dados.shape[1]
P2 = np.abs(2*Y/L)
X_ = P2[:,1:np.round(L/2)+1]
del Y, L, P2

###########################################################
#Montando vetor de caracteristicas
X = np.array([
              (dados-df)[:,0]-(dados-df)[:,1],
              (dados-df)[:,2]-(dados-df)[:,3],
              np.sum(X_[:,10:], axis=1),
              np.std((dados-df), axis=1),
              ]).T
#Separando em vetor de treino e vetor de teste
X_train, X_test, y_train, y_test = train_test_split(
         X, y, test_size=0.3, random_state=1)

#Treinando o algoritmo (ja com melhores hiperparametros)
clf = RandomForestClassifier(random_state=1, criterion='entropy',
                             max_depth=11, n_estimators=10)
clf.fit(X_train, y_train)

#Salvando Classificador Treinado
joblib.dump(clf, 'Com_Sem_Ruido.pkl')