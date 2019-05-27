# -*- coding: utf-8 -*-
import numpy as np
import math
from random import uniform
from random import randint

class Disturbios(object):
    
    def __init__(self, quant, tam, periodo, freq, snr=0):
        self.num = quant;
        self.tam = tam;
        self.periodo = periodo;
        self.freq = freq;
        self.snr = snr;
        self.t = np.linspace(0.0, self.tam/10000.0, self.tam);

    def seno_puro(self): 
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            #Valor da frequência
            freq = self.freq
            #Valor de redução na amplitude
            mult_amp = uniform(0.999, 1.001);
            seno = mult_amp*np.sin(2*np.pi*freq*self.t);
            matrix_r[i] = seno;
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
        return matrix_r;
    
    def interrupt(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            #Define onde começa e termina a interrupção
            i_begin = randint(0, self.tam-84);
            i_end = min(self.tam-1,i_begin + randint(84,self.tam));
            
            #Valor de redução na amplitude
            mult_amp = uniform(0, 0.1);
            
            #Distúrbio de interrupção
            int_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = int_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):            
            #Define onde começa e termina a interrupção
            i_begin = randint(0, self.tam-84);
            i_end = min(self.tam-1,i_begin + randint(84,self.tam));
            
            #Valor de redução na amplitude
            mult_amp = uniform(0.1, 0.9);
            
            #Distúrbio de interrupção
            sag_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = sag_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            #Define onde começa e termina a interrupção
            i_begin = randint(0, self.tam-84);
            i_end = min(self.tam-1,i_begin + randint(84,self.tam));
            
            #Valor de redução na amplitude
            mult_amp = uniform(1.1, 1.8);
            
            #Distúrbio de interrupção
            swell_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = swell_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def oscilacao(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
#            #Amplitude da oscilação (entre 0.01 e 0.3)
#            amp_osc = uniform(0.03, 0.3);
#            
#            #Parâmetros de oscilação
#            beta = randint(1, 7);
#            gama = randint(800, 900);
#            freq_o = randint(10000, 15000);
#            t1_size = 20;
#            t1 = np.linspace(0.001, 0.002, t1_size);
            
            amp_osc = uniform(0.3, 0.5)
            beta = 1
            gama = 1/uniform(0.003, 0.05)
            t1_size = randint(30,100)
            t1 = np.linspace(0.001, 0.002, t1_size)
            freq_o = randint(100, 4999)
            
            #Define onde começa e termina a interrupção
            i_begin = randint(0, self.tam-200);
            
            #Oscilação
            h1 = beta*(np.exp(-gama*t1));
            h2 = np.sin(2*math.pi*t1*freq_o);
            h3 = h1*h2;
            
            osc = np.concatenate((np.zeros(i_begin), h3, np.zeros(self.tam-(i_begin+t1_size))), axis=0);
          
            #Distúrbio de Oscilação
            osc_seno = seno + amp_osc*osc; 
            
            matrix_r[i] = osc_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def oscilacao_sem_seno(self):
        matrix_r = np.empty([self.num, self.tam])
        for i in range(0, self.num):
            
#            #Amplitude da oscilação (entre 0.01 e 0.3)
#            amp_osc = uniform(0.03, 0.3);
#            
#            #Parâmetros de oscilação
#            beta = randint(1, 7);
#            gama = randint(800, 900);
#            freq_o = randint(10000, 15000);
#            t1_size = 20;
#            t1 = np.linspace(0.001, 0.002, t1_size);
            
            amp_osc = uniform(0.3, 0.5)
            beta = 1
            gama = 1/uniform(0.003, 0.05)
            t1_size = randint(30,100)
            t1 = np.linspace(0.001, 0.002, t1_size)
            freq_o = randint(100, 4999)
            
            #Define onde começa e termina a interrupção
            i_begin = randint(0, self.tam-200);
            
            #Oscilação
            h1 = beta*(np.exp(-gama*t1));
            h2 = np.sin(2*math.pi*t1*freq_o);
            h3 = h1*h2;
            
            osc = np.concatenate((np.zeros(i_begin), h3, np.zeros(self.tam-(i_begin+t1_size))), axis=0);
          
            #Distúrbio de Oscilação
            osc_seno = amp_osc*osc; 
            
            matrix_r[i] = osc_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def flicker(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Parâmetros do flicker
            alpha = uniform(0.1, 0.2);
            #alpha = alpha if alpha != 0 else 0.01;
            freq_fck = uniform(0.1, 0.5);
            
            #Flicker
            flck = (1+(alpha*np.sin(2*np.pi*self.freq*self.t*freq_fck)));
            flck = flck*seno;
            #plt.plot(flck);
          
            #Distúrbio de Flicker            
            matrix_r[i] = flck;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def harmonica(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Número de harmônicas
            num_harm = randint(1, 6);
            harm = seno;
            
            for j in range(0, num_harm):
                #Amplitude da harmonica
                amp_h = uniform(0.05,0.2);
                #Frequência da harmônica
                freq_h = randint(2, 40);
                #Sinal com harmônica
                harm = harm + amp_h*np.sin(2*np.pi*self.t*self.freq*freq_h);
            
            #Distúrbio de Harmônicas            
            matrix_r[i] = harm;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Amplitude do distúrbio
            alpha = uniform(0.02, 0.2);
            
            #Número de distúrbios notching
            n_notch = randint(5,20);
            
            #Notching vazio
            nch = np.zeros(self.tam);
            
            for j in range(0, n_notch):
                nch[(self.tam/n_notch)*j] = 1;
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = alpha*nch + seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def notching_sem_seno(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Amplitude do distúrbio
            alpha = uniform(0.01, 0.2);
            
            #Número de distúrbios notching
            n_notch = randint(5,20);
            
            #Notching vazio
            nch = np.zeros(self.tam);
            
            for j in range(0, n_notch):
                nch[(self.tam/n_notch)*j] = 1;
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = alpha*nch;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def spike(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Amplitude do pico
            alpha = uniform(0.1, 0.5);
            
            #Onde ocorre o spike
            i_low = randint(0,480);
            i_high = min(499,i_low+randint(2,20))
            
            #Spike vazio
            spk = np.zeros(self.tam);
            spk[i_low:i_high] = 1;
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = alpha*spk + seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def spike_sem_seno(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Amplitude do pico
            alpha = uniform(0.1, 0.5);
            
            #Onde ocorre o spike
            i_low = randint(0,480);
            i_high = min(499,i_low+randint(2,20))
            
            #Spike vazio
            spk = np.zeros(self.tam);
            spk[i_low:i_high] = 1;
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = alpha*spk
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def dc_level(self):
        matrix_r = np.empty([self.num, self.tam]);
        seno = np.sin(2*np.pi*self.freq*self.t);
        for i in range(0, self.num):
            
            #Variação DC
            nivel = uniform(-0.3, 0.3);
            nivel = nivel if nivel != 0 else 0.25;
            i_begin = randint(1,50);
            i_end = i_begin + randint(150,self.tam-50);
            
            #Nivel DC
            dc = np.concatenate((seno[0:i_begin], 
                                 seno[i_begin:i_end]+nivel, 
                                 seno[i_end:self.tam]), axis=0);
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = dc;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_harm(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            harmonica = self.harmonica();
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = (sag[0] + harmonica[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_harm(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            swell = self.swell();
            harmonica = self.harmonica();
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = (swell[0] + harmonica[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def interrupt_harm(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            self.num = 1;
            seno = self.harmonica();
            seno = seno[0];
            
            #Define onde começa e termina a interrupção
            i_begin = randint(0, 300);
            i_end = min(499,i_begin + randint(84,500));
            
            #Valor de redução na amplitude
            mult_amp = uniform(0, 0.1);
            
            #Distúrbio de interrupção
            int_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = int_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def interrupt_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            self.num = 1;
            seno = self.notching();
            seno = seno[0];
            
            #Define onde começa e termina a interrupção
            i_begin = randint(0, 300);
            i_end = min(499,i_begin + randint(84,500));
            
            #Valor de redução na amplitude
            mult_amp = uniform(0, 0.1);
            
            #Distúrbio de interrupção
            int_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = int_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            notching = self.notching_sem_seno();
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = (sag[0] + notching[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            swell = self.swell();
            notching = self.notching_sem_seno();
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = (swell[0] + notching[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def interrupt_osc(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            self.num = 1;
            seno = self.oscilacao();
            seno = seno[0];
            
            #Define onde começa e termina a interrupção
            i_begin = randint(0, 300);
            i_end = min(499,i_begin + randint(84,500));
            
            #Valor de redução na amplitude
            mult_amp = uniform(0, 0.1);
            
            #Distúrbio de interrupção
            int_seno = np.concatenate((seno[0:i_begin], seno[i_begin:i_end]*mult_amp, seno[i_end:self.tam]), axis=0);  
            
            matrix_r[i] = int_seno;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_oscilacao(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            oscilacao = self.oscilacao_sem_seno();
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = (sag[0] + oscilacao[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_oscilacao(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            swell = self.swell();
            oscilacao = self.oscilacao_sem_seno();
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = (swell[0] + oscilacao[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_spike(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            spike = self.spike_sem_seno();
                    
            #Distúrbio de Spike            
            matrix_r[i] = (sag[0] + spike[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_spike(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            swell = self.swell();
            spike = self.spike_sem_seno();
                    
            #Distúrbio de Harmônicas            
            matrix_r[i] = (swell[0] + spike[0]);
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_harmonicas_oscilacao(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag + Harmonica + Oscilacao
            self.quant = 1;
            sag = self.sag();
            harmonica = self.harmonica();
            oscilacao = self.oscilacao_sem_seno();
            
                    
            #Juntando distúrbios         
            matrix_r[i] = (sag[0] + harmonica[0] + oscilacao[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_harmonicas_oscilacao(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Swell + Harmonica + Oscilacao
            self.quant = 1;
            swell = self.swell();
            harmonica = self.harmonica();
            oscilacao = self.oscilacao_sem_seno();
                    
            #Distúrbio     
            matrix_r[i] = (swell[0] + harmonica[0] + oscilacao[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def sag_harmonicas_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            sag = self.sag();
            harmonica = self.harmonica();
            notching = self.notching_sem_seno();
            
                    
            #Distúrbio            
            matrix_r[i] = (sag[0] + harmonica[0] + notching[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def swell_harmonicas_notching(self):
        matrix_r = np.empty([self.num, self.tam]);
        for i in range(0, self.num):
            
            #Sag
            self.quant = 1;
            swell = self.swell();
            harmonica = self.harmonica()
            notching = self.notching_sem_seno();
                    
            #Distúrbio        
            matrix_r[i] = (swell[0] + harmonica[0] + notching[0])/2.0;
            
            if(self.snr != 0):
                matrix_r[i] = self.noise(matrix_r[i], self.snr)
                
        return matrix_r;
    
    def noise(self, Signal, Desired_SNR_dB):
        Npts = Signal.shape[0] # Number of input time samples
        Noise = np.random.randn(Npts); # Generate initial noise; mean zero, variance one
        
        Signal_Power = np.sum(np.abs(Signal)*np.abs(Signal))/Npts;
        Noise_Power = np.sum(np.abs(Noise)*np.abs(Noise))/Npts;
        	#Initial_SNR = 10*(log10(Signal_Power/Noise_Power));
        
        K = (Signal_Power/Noise_Power)*math.pow(10,(-Desired_SNR_dB/10));  #Scale factor
        
        New_Noise = np.sqrt(K)*Noise; #Change Noise level
        	#New_Noise_Power = sum(abs(New_Noise).*abs(New_Noise))/Npts
        	#New_SNR = 10*(log10(Signal_Power/New_Noise_Power))
        
        Noisy_Signal = Signal + New_Noise;
        return Noisy_Signal