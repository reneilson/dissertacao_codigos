import numpy as np
from scipy.signal import hilbert
from scipy.ndimage import filters
import scipy.fftpack
from scipy.stats import kurtosis, skew, moment

def convert_rf(dados):
    from scipy.signal import hilbert
    from scipy.ndimage import filters
    import scipy.fftpack
    from scipy.stats import kurtosis, skew, moment

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
    df2 = filters.gaussian_filter1d(dados,10)
    
    #Calculando variaveis
    eng_data = (dados**2).sum(axis=1)
    som_max_med = np.sum(df,axis=1)*(np.max(np.abs(df), axis=1)-np.mean(abs(df), axis=1))
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
    
    X = np.array([som_max_med,
                  max_df,
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
                  std_filter_X
                  ]).T
    return X

def convert(dados):
    #Filtrando dados de entrada
    df=filters.gaussian_filter1d(dados,10)
    
    #Calculando FFT dos sinais
    Y = scipy.fftpack.fft(dados)
    L = dados.shape[1]
    P2 = np.abs(2*Y/L)
    X_ = P2[:,1:np.round(L/2)+1]
    del Y, L, P2
    
    X = np.array([
              (dados-df)[:,0]-(dados-df)[:,1],
              (dados-df)[:,2]-(dados-df)[:,3],
              np.sum(X_[:,10:], axis=1),
              np.std((dados-df), axis=1),
              ]).T
    return X

def convert_multiplos(dados):
    from scipy.ndimage import filters
    import scipy.fftpack
    #Filtrando dados de entrada
    df=filters.gaussian_filter1d(dados,10)
    
    #Calculando FFT dos sinais
    Y = scipy.fftpack.fft(dados)
    L = dados.shape[1]
    P2 = np.abs(2*Y/L)
    X_ = P2[:,1:np.round(L/2)+1]
    del Y, L, P2
    
    X = np.array([
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
    return X

def convert_sem_ruido(dados):
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
    df2 = filters.gaussian_filter1d(dados,10)
    
    #Calculando variaveis
    eng_data = (dados**2).sum(axis=1)
    som_max_med = np.sum(df,axis=1)*(np.max(np.abs(df), axis=1)-np.mean(abs(df), axis=1))
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
    
    X = np.array([#som_max_med,
              max_df,
              eng_data,
              var_X2,
              #max_mean_X2,
              #sum_X,
              #X2_sumX2,
              #max_X,
              #max_filter_X,
              sum_dados_df,
              #kurt_filter_X,
              #skew_filter_X,
              moment4_dados,
              moment2_dados,
              #skew_std_sum_X,
              #max_filter_X_sum_X,
              first_dado_df,
              #first_X_X_2,
              max_mean_X3,
              kurt_hilbert,
              #moment_2_4_hilbert,
              #moment_2_X2,
              #max_min_mean_hilbert,
              sum_filter_dados,
              std_filter_X,
              np.sum(X_[:,1:10], axis=1)*kurtosis(amplitude_envelope, axis=1),
              #np.std(X_[:,1:10], axis=1)*np.var(dados, axis=1),
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
    
    return X

def convert_dt(dados):
    Y = np.fft.fft(dados)
    L = dados.shape[1]
    P2 = np.abs(Y/L)
    X_fft = P2[:,1:np.round(L/2)+1]
    X_fft = 100*X_fft[:,10:]
    del Y, L, P2
    
    cont_mx = 0
    cont_mn = 0
    cont1 = 0
    cont2 = 0
    max_ = np.zeros(dados.shape[0])
    min_ = np.zeros(dados.shape[0])
    cont_max = np.zeros(dados.shape[0])
    cont_min = np.zeros(dados.shape[0])
    
    for i in range(X_fft.shape[0]):
        for j in range(1, X_fft.shape[1]-1):
            if ((X_fft[i, j] > X_fft[i,j+1]) and 
                (X_fft[i, j-1] < X_fft[i,j])):
                cont_mx = cont_mx + (X_fft[i,j] - X_fft[i,j-1])
                cont1 = cont1 + 1
            if ((X_fft[i, j] < X_fft[i,j+1]) and 
                (X_fft[i, j-1] > X_fft[i,j])):
                cont_mn = cont_mn + X_fft[i,j]
                cont2 = cont2 + 1
        max_[i] = cont_mx*cont1
        cont_max[i] = cont1
        min_[i] = cont_mn*cont2
        cont_min[i] = cont2
        cont_mx = 0
        cont_mn = 0
        cont1 = 0
        cont2 = 0
    
    max_aux = np.zeros(dados.shape[0])
    for i in range(dados.shape[0]):
        for j in range(1, dados.shape[1]-1):
            if ((dados[i, j] > dados[i,j+1]) and 
                (dados[i, j-1] < dados[i,j])):
                max_aux[i] = max_aux[i] + 1
    max_aux = max_aux-3
    
    del i, j, cont1
    del cont2, cont_mn, cont_mx
    
    aux = np.zeros(dados.shape)
    for i in range(dados.shape[0]):
        for j in range(dados.shape[1]):
            aux[i,j] = np.max(abs(dados[i,j:j+100%dados.shape[1]]))
            
    aux2 = np.zeros(X_fft.shape)
    for i in range(X_fft.shape[0]):
        for j in range(X_fft.shape[1]):
            aux2[i,j] = np.max(abs(X_fft[i,j:j+50%X_fft.shape[1]]))
            
    tam = dados.shape[1]
    #std = np.std(np.array([np.mean(dados[:,0:tam/5], axis=1), 
    #              np.mean(dados[:,tam/5:(tam/5)*2], axis=1),
    #              np.mean(dados[:,(tam/5)*2:(tam/5)*3], axis=1), 
    #              np.mean(dados[:,(tam/5)*3:(tam/5)*4], axis=1), 
    #              np.mean(dados[:,(tam/5)*4:], axis=1)]).T,axis=1)
            
    
    maxi = 4
    std = np.zeros([maxi, dados.shape[0]])
    for i in range(maxi):
        std[i,:] = np.mean(dados[:,i*(tam/maxi):(i+1)*(tam/maxi)],axis=1)
    
    soma_fft = np.sum(X_fft, axis=1)
    rms_dados = np.sqrt(np.mean(dados**2, axis=1))
    dist_max_mean_dados = np.max(np.abs(dados), axis=1)-np.mean(dados, axis=1)
    dist_min_mean_dados = -np.min(dados, axis=1)+np.mean(dados, axis=1)
    dist_cont = cont_min-cont_max
    min_env_dados = np.min(aux[:,:450], axis=1)
    max_env_dados = np.max(aux[:,:450], axis=1)
    dist_env_fft_mx = np.max(aux2, axis=1)-np.min(aux2[:,:200], axis=1)
    dist_env_fft_mx_mean = np.max(aux2, axis=1)-np.mean(aux2, axis=1)
    std_ = np.std(std.T,axis=1)
    
    X = np.array([soma_fft,
                  max_**2,  
                  rms_dados,
                  std_,
                  dist_max_mean_dados,
                  dist_min_mean_dados,
                  dist_cont,
                  min_env_dados,
                  max_env_dados,
                  dist_env_fft_mx,
                  dist_env_fft_mx_mean,
                  max_aux
                  ]).T
    return X

def convert_com_ruido(dados):
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
    som_max_med = np.sum(df,axis=1)*(np.max(np.abs(df), axis=1)-np.mean(abs(df), axis=1))
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
    
    X = np.array([#som_max_med,
              max_df,
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
    
    return X

def convert_new(dados):
    #Filtrando dados de entrada
    df=filters.gaussian_filter1d(dados,5)
    
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
    #var_X2 = np.var(X_2[:,2:10],axis=1)
    #max_mean_X2 = np.max(X_[:,3:], axis=1)-np.mean(X_[:,3:], axis=1)
    #sum_X = np.sum((100*X_[:,3:])**4, axis=1)
    #X2_sumX2 = (X_2[:,2])*np.sum(X_2[:,:10], axis=1)
    #max_X = np.max(X_[:,3:], axis=1)
    #max_filter_X = np.max(df2[:,50:],axis=1)
    sum_dados_df = np.sum((dados-df)**2, axis=1)
    kurt_filter_X = kurtosis(df2, axis=1)
    skew_filter_X = skew(df2, axis=1)
    moment4_dados = moment(dados, 4, axis=1)
    moment2_dados = moment(dados, 2, axis=1)
    #skew_std_sum_X = skew(X_[:,4:], axis=1)+1000*np.std(X_[:,3:],axis=1)+10*np.sum(X_[:,4:], axis=1)
    #max_filter_X_sum_X = np.max(df2[:,50:],axis=1)*np.sum((100*X_[:,3:])**4, axis=1)
    #first_dado_df = (dados-df)[:,0]
    first_X_X_2 = X_[:,2]+X_2[:,2]
    max_mean_X3 = np.max(X_3[:,4:], axis=1)-np.mean(X_3, axis=1)
    kurt_hilbert = kurtosis(amplitude_envelope, axis=1)
    moment_2_4_hilbert = moment(amplitude_envelope,  2, axis=1)*moment(amplitude_envelope,  4, axis=1)
    #moment_2_X2 = moment(X_2[:,3:], 2, axis=1)
    #max_min_mean_hilbert = (np.max(amplitude_envelope, axis=1) +np.min(amplitude_envelope, axis=1))*np.mean(amplitude_envelope, axis=1)
    #sum_filter_dados = np.sum(df, axis=1)
    #std_filter_X = np.std(df2, axis=1)
    
    
    aux = np.zeros(dados.shape)
    for i in range(dados.shape[0]):
        for j in range(dados.shape[1]):
            aux[i,j] = np.max(abs(dados[i,j:j+100%dados.shape[1]]))
            
    min_env_dados = np.min(aux[:,:450], axis=1)
    max_env_dados = np.max(aux[:,:450], axis=1)
    
    ############################################################
    #Montando Vetor de Caracteristicas
    X = np.array([max_df,
                  eng_data,
                  #np.max(amplitude_envelope, axis=1),
                  #np.min(amplitude_envelope, axis=1),
                  #var_X2,
                  min_env_dados,
                  max_env_dados,
                  #max_mean_X2,
                  #sum_X,
                  #X2_sumX2,
                  #max_X,
                  #max_filter_X,
                  sum_dados_df,
                  kurt_filter_X,
                  skew_filter_X,
                  moment4_dados,
                  moment2_dados,
                  #skew_std_sum_X,
                  #max_filter_X_sum_X,
                  #first_dado_df,
                  first_X_X_2,
                  max_mean_X3,
                  kurt_hilbert,
                  moment_2_4_hilbert,
                  #moment_2_X2,
                  #max_min_mean_hilbert,
                  #sum_filter_dados,
                  #std_filter_X,
                  #np.sum(X_[:,1:10], axis=1)*kurtosis(amplitude_envelope, axis=1),
                  #np.std(X_[:,1:10], axis=1)*np.var(dados, axis=1),
                  #np.sum(X_[:,10:], axis=1),
                  #np.std(X_[:,10:], axis=1),
                  np.max(X_[:,10:], axis=1)-np.mean(X_[:,10:], axis=1),
                  np.max(X_[:,50:100], axis=1)-np.mean(X_[:,50:100], axis=1),
                  np.max(X_[:,150:200], axis=1)-np.mean(X_[:,150:200], axis=1),
                  np.max(X_[:,200:], axis=1)-np.mean(X_[:,200:], axis=1),
                  #np.sum(X_[:,:50], axis=1)-
                  #np.sum(X_[:,50:100], axis=1)*np.sum(X_[:,150:200], axis=1),
                  #np.sum(X_[:,200:], axis=1),
                  np.std(X_[:,:50], axis=1),
                  np.std(X_[:,150:200], axis=1),
                  #np.sum(X_[:,0:3],axis=1),
                  dados[:,0],
                  df[:,0],
                  #np.sum(df2, axis=1)*skew(X_, axis=1),
                  np.std(dados-df, axis=1),
                  #(kurtosis(X_, axis=1) - kurtosis(df2, axis=1))*kurtosis(X_-df2, axis=1),
                  np.sqrt(np.mean(dados**2, axis=1)),
                  np.sqrt(np.mean(df**2, axis=1)),
                  np.sqrt(np.mean(amplitude_envelope**2, axis=1)),
                  np.std(amplitude_envelope, axis=1),
                  np.sum(amplitude_envelope, axis=1),
                  #np.max(amplitude_envelope, axis=1),
                  ]).T
    return X