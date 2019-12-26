import numpy as np
from scipy.signal import hilbert
from scipy.ndimage import filters
import scipy.fftpack
from scipy.stats import kurtosis, skew, moment

def convert(dados):
    #Data filtering
    df=filters.gaussian_filter1d(dados,5)
    
    #FFT Calculation
    Y = scipy.fftpack.fft(dados)
    L = dados.shape[1]
    P2 = np.abs(2*Y/L)
    X_ = P2[:,1:np.round(L/2)+1]
    del Y, L, P2
    
    #FFT from filtered signals
    Y = np.fft.fft(df)
    L = dados.shape[1]
    P2 = np.abs(2*Y/L)
    X_2 = P2[:,1:np.round(L/2)+1]
    del Y, L, P2
    
    #FFT from diference of filtered and input signals
    Y = np.fft.fft(dados-df)
    L = (dados-df).shape[1]
    P2 = np.abs(Y/L)
    X_3 = P2[:,1:np.round(L/2)+1]
    del Y, L, P2
    
    #Hilbert Transform
    analytic_signal = hilbert(dados)
    amplitude_envelope = np.abs(analytic_signal)
    
    #Filtering FFT
    df2 = filters.gaussian_filter1d(X_,10)
    
    #Applying some statistical calculations
    eng_data = (dados**2).sum(axis=1)
    max_df = np.max(abs(df), axis=1)
    sum_dados_df = np.sum((dados-df)**2, axis=1)
    kurt_filter_X = kurtosis(df2, axis=1)
    skew_filter_X = skew(df2, axis=1)
    moment4_dados = moment(dados, 4, axis=1)
    moment2_dados = moment(dados, 2, axis=1)
    first_X_X_2 = X_[:,2]+X_2[:,2]
    max_mean_X3 = np.max(X_3[:,4:], axis=1)-np.mean(X_3, axis=1)
    kurt_hilbert = kurtosis(amplitude_envelope, axis=1)
    moment_2_4_hilbert = moment(amplitude_envelope,  2, axis=1)*moment(amplitude_envelope,  4, axis=1)
    
    #Mathematical Morfology from Input signal
    aux = np.zeros(dados.shape)
    for i in range(dados.shape[0]):
        for j in range(dados.shape[1]):
            aux[i,j] = np.max(abs(dados[i,j:j+100%dados.shape[1]]))        
    min_env_dados = np.min(aux[:,:450], axis=1)
    max_env_dados = np.max(aux[:,:450], axis=1)
    
    ############################################################
    #Final Feature Vector
    X = np.array([max_df,
                  eng_data,
                  min_env_dados,
                  max_env_dados,
                  sum_dados_df,
                  kurt_filter_X,
                  skew_filter_X,
                  moment4_dados,
                  moment2_dados,
                  first_X_X_2,
                  max_mean_X3,
                  kurt_hilbert,
                  moment_2_4_hilbert,
                  np.max(X_[:,10:], axis=1)-np.mean(X_[:,10:], axis=1),
                  np.max(X_[:,50:100], axis=1)-np.mean(X_[:,50:100], axis=1),
                  np.max(X_[:,150:200], axis=1)-np.mean(X_[:,150:200], axis=1),
                  np.max(X_[:,200:], axis=1)-np.mean(X_[:,200:], axis=1),
                  np.std(X_[:,:50], axis=1),
                  np.std(X_[:,150:200], axis=1),
                  dados[:,0],
                  df[:,0],
                  np.std(dados-df, axis=1),
                  np.sqrt(np.mean(dados**2, axis=1)),
                  np.sqrt(np.mean(df**2, axis=1)),
                  np.sqrt(np.mean(amplitude_envelope**2, axis=1)),
                  np.std(amplitude_envelope, axis=1),
                  np.sum(amplitude_envelope, axis=1),
                  ]).T
    return X
