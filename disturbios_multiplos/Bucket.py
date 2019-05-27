# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 12:29:40 2017

@author: Rene
"""

from disturbios import Disturbios
import numpy as np
from random import uniform

class Bucket(object):
    
    def __init__(self, n_ex, snr, freq=60, simples=True):
        self.n_ex = n_ex;
        self.snr = snr;
        self.freq = freq 
        self.simples = simples;
        
    def createBucket(self):
        dist = []
        for i in range(self.n_ex):
            freq = self.freq + uniform(-0.01*self.freq, 0.01*self.freq);
            d = Disturbios(1, 500, 3, freq, self.snr);
            
            if(self.simples):
                #Distúrbios Simples
                sn = d.seno_puro();
                i = d.interrupt();
                sg = d.sag();
                sw = d.swell();
                osc = d.oscilacao();
                flk = d.flicker();
                hm = d.harmonica();
                ntc = d.notching();
                spk = d.spike();
                dcl = d.dc_level();
                
                if dist == []:
                    dist = np.concatenate((sn, i, sg, sw, osc, flk, hm, ntc, 
                                           spk, dcl
                                           ), axis=0);
                else:
                    dist = np.concatenate((dist, sn, i, sg, sw, osc, flk, hm, 
                                           ntc, spk, dcl, 
                                           ), axis=0);
            else:
                #Distúrbios Simples
                sn = d.seno_puro();
                i = d.interrupt();
                sg = d.sag();
                sw = d.swell();
                osc = d.oscilacao();
                flk = d.flicker();
                hm = d.harmonica();
                ntc = d.notching();
                spk = d.spike();
                dcl = d.dc_level();
                #Distúrbios Multiplos
                swh = d.swell_harm();
                sgh = d.sag_harm();
                swn = d.swell_notching();
                sgn = d.sag_notching();
                swo = d.swell_oscilacao();
                sgo = d.sag_oscilacao();
                sgs = d.sag_spike();
                sws = d.swell_spike();
                sghn = d.sag_harmonicas_notching();
                swhn = d.swell_harmonicas_notching();
                sgho = d.sag_harmonicas_oscilacao();
                swho = d.swell_harmonicas_oscilacao();

                if dist == []:
                    dist = np.concatenate((sn, i, sg, sw, osc, flk, hm, ntc, 
                                           spk, dcl, swh, sgh, swn, sgn, swo, 
                                           sgo, sgs, sws, sghn, swhn, sgho, swho), 
                                           axis=0);
                else:
                    dist = np.concatenate((dist, sn, i, sg, sw, osc, flk, hm, ntc, 
                                           spk, dcl, swh, sgh, swn, sgn, swo, 
                                           sgo, sgs, sws, sghn, swhn, sgho, swho), 
                                           axis=0);
            #END IF-ELSE
        return dist;
    
    def createLabels(self):
        #Criando label (n_d representa o número de distúrbios)
        if(self.simples):
            n_d = 10;
        else:
            n_d = 22; 
        lbl = np.zeros([self.n_ex*n_d, n_d]);
        
        for i in range(0,(self.n_ex*n_d)):
            lbl[i][i//self.n_ex] = 1;
        
        return lbl;
    
    def createLabelsSVM(self):
        #Criando label (n_d representa o número de distúrbios)
        if(self.simples):
            n_d = 10;
        else:
            n_d = 22; 
        lbl = np.zeros([self.n_ex*n_d]);
        
        for i in range(0,(self.n_ex*n_d)):
            lbl[i]=i%n_d;
            
        return lbl;
    
    def createFile(self, fileNameBucket, fileNameLabel, fileNameLabelSVM):
        np.save(fileNameBucket, self.createBucket());
        np.save(fileNameLabel, self.createLabels());
        np.save(fileNameLabelSVM, self.createLabelsSVM());
        
        