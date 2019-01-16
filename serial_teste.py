import serial
import numpy as np

porta = '/dev/ttyUSB0'
baud_rate = 115200

for k in range(30):
    #Ler base de dados com 200 disturbios simples
    d = np.load('dados_'+str(k)+'.npy')
    #Variavel para capturar classes
    #cl = np.zeros(200)
    #Variavel para capturar tempo de processamento
    tm = np.zeros(200)
    #Sao mandados de um em um disturbio presente na base
    for j in range(0,200):
        #abre a porta serial
        Obj_porta = serial.Serial(porta, baud_rate)
        #envia amostra por amostra
        for i in range(512):
            #insere o \n ao final
            valor = str(d[j,i]) + b"\n"
            try:
                #tenta enviar o dado pela porta serial
                Obj_porta.write(valor)
                #recebe o valor enviado (serve como garantia de envio)
                valor = float(Obj_porta.readline().splitlines()[0])
            except serial.SerialException:
                print "ERRO" 
        #apos todos os dados enviados, aguarda o resultado (classe referente)
        valor = Obj_porta.readline()
        #print "CLASSE: ", valor
        #cl[j] = int(valor.splitlines()[0])
        #recebe o tempo de processamento via porta serial
        valor = Obj_porta.readline()
        #print "TIME:", valor
        #armazena tempo de processamento
        tm[j] = int(valor.splitlines()[0])
        #fecha a porta serial
        Obj_porta.close()
    #apos todos os 200 disturbios classificados, armazena os tempos e arquivo
    np.save('time' + str(k), tm)
    #np.save('cl' + str(k), cl)
    #print("END " + str(k))