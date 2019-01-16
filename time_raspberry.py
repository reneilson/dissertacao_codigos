import numpy as np
from Bucket import Bucket
from classifier import classifier
import time


for ii in range(30):    
    #leitura dos dados
    dados = np.load('dados' +str(ii) + '.npy')
    
    #inicia o tempo utilizado para realizar a classificacao
    inicio = time.time()
    
    #Calculo da FFT
    Y = np.fft.fft(dados)
    L = dados.shape[1]
    P2 = np.abs(Y/L)
    X_fft = P2[:,1:np.round(L/2)+1]
    X_fft = 100*X_fft[:,10:]
    del Y, L, P2
    eng_data = (dados**2).sum(axis=1)
    
    cont_mx = 0
    cont_mn = 0
    cont1 = 0
    cont2 = 0
    max_ = np.zeros(dados.shape[0])
    min_ = np.zeros(dados.shape[0])
    cont_max = np.zeros(dados.shape[0])
    cont_min = np.zeros(dados.shape[0])
    cont_max_dados = np.zeros(dados.shape[0])
    
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
    maxi = 4
    std = np.zeros([maxi, dados.shape[0]])
    for i in range(maxi):
        std[i,:] = np.mean(dados[:,i*(tam/maxi):(i+1)*(tam/maxi)],axis=1)
    
    #caracteristicas utilizadas
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
    
    cl = np.zeros([200])
    
    #classificacao dos 200 dados
    for i in range(dados.shape[0]):
        cl[i] = classifier(soma_fft[i], max_[i], rms_dados[i], 
                   std_[i], dist_max_mean_dados[i],
                   dist_min_mean_dados[i], dist_cont[i],
                   min_env_dados[i], max_env_dados[i],
                   dist_env_fft_mx[i], dist_env_fft_mx_mean[i],
                   max_aux[i])
        
    #finaliza o tempo de processamento
    fim = time.time()
    #np.save('cl_rasp' + str(ii), cl)
    #Salva o tempo final 
    np.save('time' + str(ii), fim-inicio)
    
