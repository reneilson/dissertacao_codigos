#include <math.h>
#include "arduinoFFT.h"
#define SAMPLING_FREQUENCY 10000 //Hz, must be less than 10000 due to ADC
 
arduinoFFT FFT = arduinoFFT();

unsigned int sampling_period_us;
unsigned long microseconds;
 
const int SOMA = 0;
const int MEDIA = 1;
const int ENERGIA = 2;
const int NUM_MAX = 3;
const int VAL_MAX = 4;
const int NUM_MIN = 5;
const int VAL_MIN = 6;
const int RMS = 7;
const int MAX_ = 8;
const int MIN_ = 9;
const int SAMPLES = 512;
const int SAMPLES_DATA = 512;

void func_vector_fft(float* vector, int tam, double *dados)
{
  int i;

  //Vetor de retorno: 
  //0- soma, 1- media, 2- energia
  //3- numero de maximos, 4- valores maximos
  //5- numero de minimos, 6- valores minimos
  //7- rms
  //8- maximo
  //9- minimo
  
  for (i = 0; i < 10; i++){
    vector[i] = 0;
  }
  
  for (i = 1; i < tam; i++){
    vector[0] = vector[0] + dados[i];  
    vector[2] = vector[2] + dados[i]*dados[i];
    if ((dados[i] > dados[i+1]) && (dados[i] > dados[i-1])){
      vector[3]++;
      vector[4] = vector[4] + (dados[i] - dados[i-1]);
    }
    if ((dados[i] < dados[i+1]) && (dados[i] < dados[i-1])){
      vector[5]++;
      vector[6] = vector[6] + dados[i];
    }
    if (dados[i] > vector[8])
      vector[8] = dados[i];
    if (dados[i] < vector[9])
      vector[9] = dados[i];
  }

  vector[1] = vector[0]/tam;
  vector[7] = sqrt(vector[2]/tam);
  
}

void func_vector(float* vector, int tam, float *dados)
{
  int i;

  //Vetor de retorno: 
  //0- soma, 1- media, 2- energia
  //3- numero de maximos, 4- valores maximos
  //5- numero de minimos, 6- valores minimos
  //7- rms
  //8- maximo
  //9- minimo
  
  for (i = 0; i < 10; i++){
    vector[i] = 0;
  }
  
  for (i = 0; i < tam; i++){
    vector[0] = vector[0] + dados[i];  
    vector[2] = vector[2] + dados[i]*dados[i];
    if ((i > 0 && dados[i] > dados[i+1]) && (dados[i] > dados[i-1])){
      vector[3]++;
      vector[4] = vector[4] + (dados[i] - dados[i-1]);
    }
    if ((i > 0 && dados[i] < dados[i+1]) && (dados[i] < dados[i-1])){
      vector[5]++;
      vector[6] = vector[6] + dados[i];
    }
    if (dados[i] > vector[8])
      vector[8] = dados[i];
    if (dados[i] < vector[9])
      vector[9] = dados[i];
  }

  vector[1] = vector[0]/tam;
  vector[7] = sqrt(vector[2]/tam);
  
}

float std_envelop(int tam, float *dados){
  float vector[]= {0,0,0,0};
  int i, j;
  float soma = 0, media =0, std = 0;

  for (j = 0; j < 4; j++){
    for (i = j*(tam/4); i < (j+1)*(tam/4); i++){
      soma = soma + dados[i];
    }
    vector[j] = soma / (tam/4);
    soma = 0;
    media = media + vector[j];
  } 
  media = media/4;

  for (i = 0; i < 4; i++){
    std = std + pow(vector[i] - media,2);
  } 

  std = sqrt(std/4);
  return std;
}

void extremos(float* max_min, int mn, int mx, int tam, float* dados){
  int i;
  max_min[0] = 0;
  max_min[1] = 100;
  max_min[2] = 0;
  for(i = mn; i < mx; i++){
    if(i >= tam)
      break;
    max_min[0] = max(fabs(dados[i]),fabs(max_min[0]));
    max_min[1] = min(max_min[1],dados[i]);
    max_min[2] = max_min[2] + dados[i];
  }
  max_min[2] = max_min[2]/(mx-mn);
}

void envelop(float* v_envelop, int tam, int w, float* dados){
  int i;
  float* max_min = (float *)malloc(sizeof(float) * 3);
  for (i = 0; i < tam; i++){
    extremos(max_min, i, i+w, tam, dados);
    v_envelop[i] = max_min[0];
  }
}

void extremos_fft(float* max_min, int mn, int mx, int tam, double* dados){
  max_min[0] = 0;
  max_min[1] = 100;
  max_min[2] = 0;
  for(int i = mn; i < mx; i++){
    if(i >= tam)
      break;
    max_min[0] = max(fabs(max_min[0]), fabs(dados[i]));
    if (dados[i] < max_min[1])
      max_min[1] = dados[i];
    max_min[2] = max_min[2] + dados[i];
  }
  max_min[2] = max_min[2]/(mx-mn);
}

void envelop_fft(float* v_envelop, int tam, int w, double* dados){
  int i;
  float* max_min = (float *)malloc(sizeof(float) * 3);
  for (i = 0; i < tam; i++){
    extremos_fft(max_min, i, i+w, tam, dados);
    v_envelop[i] = max_min[0];
  }
}

void construct_arg(float* X, float* dados){
  
  for(int i=0; i<12; i++){
    X[i] = 0.0;
  }  
  double* X_fft = (double *)malloc(sizeof(double) * SAMPLES);
  double* vImag = (double *)malloc(sizeof(double) * SAMPLES);

  for(int i = 0; i < SAMPLES_DATA; i++){
    X_fft[i] = dados[i];
    vImag[i] = 0;
  }
  
  //FFT.Windowing(X_fft, SAMPLES, FFT_WIN_TYP_HAMMING, FFT_FORWARD);
  FFT.Compute(X_fft, vImag, SAMPLES, FFT_FORWARD);
  FFT.ComplexToMagnitude(X_fft, vImag, SAMPLES);

  free(vImag);
  vImag = NULL;
  
  double* _X_fft = (double *)malloc(sizeof(double) * (SAMPLES - 10));
  
  for(int i =11; i < (SAMPLES/2)+1; i++){
    _X_fft[i-11] = X_fft[i]/(SAMPLES);
    _X_fft[i-11] *= 100;
  }

  free(X_fft);
  X_fft = NULL;
  
  int tam_fft = (SAMPLES/2)-10;
  
  float* vector_fft = (float *)malloc(sizeof(float) * 10);
  func_vector_fft(vector_fft, tam_fft, _X_fft);
  
  X[0] = vector_fft[SOMA];
  X[1] = vector_fft[VAL_MAX]*vector_fft[NUM_MAX];
  X[6] = vector_fft[NUM_MAX] - vector_fft[NUM_MIN];
  
  free(vector_fft);
  vector_fft = NULL;
  
  float* envelop_fft_data = (float *)malloc(sizeof(float) * tam_fft);
  envelop_fft(envelop_fft_data, tam_fft, 50, _X_fft);
  float* max_min_fft = (float *)malloc(sizeof(float) * 3);
  extremos(max_min_fft, 0, tam_fft, tam_fft, envelop_fft_data);
  
  X[9] = max_min_fft[0];
  X[10] = max_min_fft[0] - max_min_fft[2];
  
  extremos(max_min_fft, 0, 200, tam_fft, envelop_fft_data);
  X[9] = X[9] - max_min_fft[1];

  free(envelop_fft_data);
  free(max_min_fft);
  free(_X_fft);
  max_min_fft = NULL;
  envelop_fft_data = NULL;
  _X_fft = NULL;
  
  float* vector_dados = (float *)malloc(sizeof(float) * 10);
  func_vector(vector_dados, SAMPLES, dados);
  
  X[2] = vector_dados[RMS];
  if(vector_dados[MAX_] < -vector_dados[MIN_])
    vector_dados[MAX_] = -vector_dados[MIN_];
  X[4] = vector_dados[MAX_] - vector_dados[MEDIA];
  X[5] = vector_dados[MEDIA] - vector_dados[MIN_];
  X[11] = vector_dados[NUM_MAX]-3;

  free(vector_dados);
  vector_dados = NULL;

  X[3] = std_envelop(SAMPLES, dados);
  
  float* envelop_dados = (float *)malloc(sizeof(float) * SAMPLES);
  envelop(envelop_dados, SAMPLES, 100, dados);
  float* max_min_dados = (float *)malloc(sizeof(float) * 3);
  extremos(max_min_dados, 0, 450, SAMPLES, envelop_dados);
  
  X[7] = max_min_dados[1];
  X[8] = max_min_dados[0];

  free(envelop_dados);
  free(max_min_dados);
  max_min_dados = NULL;
  envelop_dados = NULL;
}

int classifier(float soma_fft, float max_, float rms_dados, 
               float std_, float dist_max_mean_dados,
               float dist_min_mean_dados, float dist_cont,
               float min_env_dados, float max_env_dados,
               float dist_env_fft_mx, float dist_env_fft_mx_mean,
               float max_aux){
    if ( min_env_dados <= 0.9937914609909058 ){
        if ( rms_dados <= 0.6837615966796875 ){
            if ( min_env_dados <= 0.18013443052768707 ){
                if ( min_env_dados <= 0.10111333429813385 ){
                    return 1;
                }else{
                    if ( max_aux <= 33.5 ){
                        if ( dist_min_mean_dados <= 1.095322847366333 ){
                            if ( min_env_dados <= 0.10353253036737442 ){
                                if ( max_env_dados <= 1.0010888576507568 ){
                                    return 2;
                                }else{
                                    return 1;
                                }
                            }else{
                                return 2;
                            }
                        }else{
                            return 1;
                        }
                    }else{
                        if ( max_env_dados <= 1.0706720352172852 ){
                            if ( min_env_dados <= 0.14462760090827942 ){
                                if ( max_env_dados <= 1.010812520980835 ){
                                    if ( std_ <= 0.15435299277305603 ){
                                        return 1;
                                    }else{
                                        return 2;
                                    }
                                }else{
                                    if ( dist_max_mean_dados <= 1.039196252822876 ){
                                        return 1;
                                    }else{
                                        if ( std_ <= 0.1986302137374878 ){
                                            return 1;
                                        }else{
                                            if ( soma_fft <= 52.92010498046875 ){
                                                return 2;
                                            }else{
                                                return 1;
                                            }
                                        }
                                    }
                                }
                            }else{
                                if ( max_ <= 1183545.75 ){
                                    if ( max_env_dados <= 1.0101730823516846 ){
                                        return 1;
                                    }else{
                                        if ( max_aux <= 153.0 ){
                                            return 2;
                                        }else{
                                            if ( min_env_dados <= 0.15602317452430725 ){
                                                return 1;
                                            }else{
                                                return 2;
                                            }
                                        }
                                    }
                                }else{
                                    return 1;
                                }
                            }
                        }else{
                            return 1;
                        }
                    }
                }
            }else{
                if ( rms_dados <= 0.6528249382972717 ){
                    if ( dist_min_mean_dados <= 1.0691787004470825 ){
                        if ( dist_max_mean_dados <= 1.0610541105270386 ){
                            if ( soma_fft <= 58.99302673339844 ){
                                return 2;
                            }else{
                                if ( max_ <= 583700.25 ){
                                    return 1;
                                }else{
                                    if ( soma_fft <= 69.84739685058594 ){
                                        return 2;
                                    }else{
                                        return 1;
                                    }
                                }
                            }
                        }else{
                            if ( rms_dados <= 0.6299293041229248 ){
                                if ( min_env_dados <= 0.24093642830848694 ){
                                    if ( max_aux <= 159.0 ){
                                        if ( dist_env_fft_mx_mean <= 0.602356493473053 ){
                                            return 1;
                                        }else{
                                            return 2;
                                        }
                                    }else{
                                        return 1;
                                    }
                                }else{
                                    return 2;
                                }
                            }else{
                                if ( rms_dados <= 0.6356284618377686 ){
                                    return 1;
                                }else{
                                    if ( max_aux <= 17.0 ){
                                        return 2;
                                    }else{
                                        if ( dist_env_fft_mx_mean <= 1.0999882221221924 ){
                                            if ( max_ <= 1377758.75 ){
                                                return 2;
                                            }else{
                                                return 1;
                                            }
                                        }else{
                                            return 1;
                                        }
                                    }
                                }
                            }
                        }
                    }else{
                        if ( min_env_dados <= 0.28379958868026733 ){
                            if ( rms_dados <= 0.5247938632965088 ){
                                if ( dist_env_fft_mx <= 0.28757375478744507 ){
                                    return 1;
                                }else{
                                    return 2;
                                }
                            }else{
                                if ( dist_min_mean_dados <= 1.1720980405807495 ){
                                    if ( std_ <= 0.2353316694498062 ){
                                        if ( dist_env_fft_mx <= 1.6502805948257446 ){
                                            if ( dist_max_mean_dados <= 0.9575475454330444 ){
                                                return 1;
                                            }else{
                                                return 2;
                                            }
                                        }else{
                                            if ( dist_min_mean_dados <= 1.1212985515594482 ){
                                                return 1;
                                            }else{
                                                return 1;
                                            }
                                        }
                                    }else{
                                        return 1;
                                    }
                                }else{
                                    return 1;
                                }
                            }
                        }else{
                            if ( rms_dados <= 0.6327909231185913 ){
                                return 2;
                            }else{
                                if ( dist_max_mean_dados <= 0.9144288301467896 ){
                                    return 1;
                                }else{
                                    if ( dist_env_fft_mx <= 1.6304149627685547 ){
                                        if ( min_env_dados <= 0.38052141666412354 ){
                                            if ( dist_min_mean_dados <= 1.1173427104949951 ){
                                                return 2;
                                            }else{
                                                return 1;
                                            }
                                        }else{
                                            return 2;
                                        }
                                    }else{
                                        return 1;
                                    }
                                }
                            }
                        }
                    }
                }else{
                    if ( dist_env_fft_mx <= 1.0390160083770752 ){
                        if ( soma_fft <= 31.34770965576172 ){
                            return 2;
                        }else{
                            if ( soma_fft <= 72.73396301269531 ){
                                if ( dist_min_mean_dados <= 1.0934584140777588 ){
                                    if ( soma_fft <= 35.427452087402344 ){
                                        if ( rms_dados <= 0.6727970838546753 ){
                                            if ( min_env_dados <= 0.9110903143882751 ){
                                                return 2;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            if ( rms_dados <= 0.6751012802124023 ){
                                                return 4;
                                            }else{
                                                return 2;
                                            }
                                        }
                                    }else{
                                        if ( max_ <= 1649208.25 ){
                                            if ( dist_min_mean_dados <= 1.0106229782104492 ){
                                                return 2;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            return 2;
                                        }
                                    }
                                }else{
                                    return 2;
                                }
                            }else{
                                if ( min_env_dados <= 0.958419919013977 ){
                                    return 4;
                                }else{
                                    return 2;
                                }
                            }
                        }
                    }else{
                        if ( rms_dados <= 0.6792887449264526 ){
                            if ( soma_fft <= 35.844749450683594 ){
                                if ( max_aux <= 6.5 ){
                                    return 2;
                                }else{
                                    return 4;
                                }
                            }else{
                                if ( min_env_dados <= 0.8622983694076538 ){
                                    return 4;
                                }else{
                                    if ( dist_min_mean_dados <= 1.0079259872436523 ){
                                        return 4;
                                    }else{
                                        return 2;
                                    }
                                }
                            }
                        }else{
                            if ( dist_env_fft_mx <= 1.3367085456848145 ){
                                return 4;
                            }else{
                                if ( dist_min_mean_dados <= 1.0159471035003662 ){
                                    return 4;
                                }else{
                                    return 8;
                                }
                            }
                        }
                    }
                }
            }
        }else{
            if ( dist_min_mean_dados <= 1.0916576385498047 ){
                if ( soma_fft <= 46.16716003417969 ){
                    if ( max_env_dados <= 1.0071170330047607 ){
                        if ( dist_min_mean_dados <= 1.0124049186706543 ){
                            if ( dist_max_mean_dados <= 0.9815332889556885 ){
                                return 6;
                            }else{
                                return 4;
                            }
                        }else{
                            if ( dist_env_fft_mx_mean <= 0.7235246300697327 ){
                                return 2;
                            }else{
                                if ( rms_dados <= 0.6884053945541382 ){
                                    if ( dist_env_fft_mx <= 1.2838664054870605 ){
                                        return 4;
                                    }else{
                                        if ( std_ <= 0.21445265412330627 ){
                                            return 8;
                                        }else{
                                            return 4;
                                        }
                                    }
                                }else{
                                    return 8;
                                }
                            }
                        }
                    }else{
                        if ( soma_fft <= 36.83250427246094 ){
                            if ( rms_dados <= 0.6965171098709106 ){
                                if ( soma_fft <= 33.497764587402344 ){
                                    if ( rms_dados <= 0.6863300800323486 ){
                                        if ( max_env_dados <= 1.099808692932129 ){
                                            return 2;
                                        }else{
                                            return 9;
                                        }
                                    }else{
                                        return 9;
                                    }
                                }else{
                                    if ( std_ <= 0.1865592896938324 ){
                                        return 9;
                                    }else{
                                        if ( dist_max_mean_dados <= 1.0398244857788086 ){
                                            if ( dist_env_fft_mx <= 0.8315548896789551 ){
                                                return 8;
                                            }else{
                                                return 9;
                                            }
                                        }else{
                                            return 4;
                                        }
                                    }
                                }
                            }else{
                                if ( max_ <= 0.00017077763914130628 ){
                                    return 5;
                                }else{
                                    if ( soma_fft <= 29.143638610839844 ){
                                        return 5;
                                    }else{
                                        return 9;
                                    }
                                }
                            }
                        }else{
                            if ( min_env_dados <= 0.8571571111679077 ){
                                if ( rms_dados <= 0.6914154291152954 ){
                                    return 4;
                                }else{
                                    return 9;
                                }
                            }else{
                                if ( max_env_dados <= 1.048864722251892 ){
                                    return 8;
                                }else{
                                    return 4;
                                }
                            }
                        }
                    }
                }else{
                    if ( std_ <= 0.19409391283988953 ){
                        if ( dist_env_fft_mx_mean <= 1.5470330715179443 ){
                            return 0;
                        }else{
                            return 8;
                        }
                    }else{
                        if ( dist_env_fft_mx_mean <= 1.6095737218856812 ){
                            if ( dist_cont <= 0.5 ){
                                return 9;
                            }else{
                                return 8;
                            }
                        }else{
                            return 6;
                        }
                    }
                }
            }else{
                if ( std_ <= 0.21673263609409332 ){
                    if ( std_ <= 0.19221431016921997 ){
                        if ( dist_max_mean_dados <= 1.176797866821289 ){
                            if ( dist_env_fft_mx <= 0.7322201728820801 ){
                                if ( dist_env_fft_mx <= 0.36461615562438965 ){
                                    return 2;
                                }else{
                                    return 9;
                                }
                            }else{
                                if ( soma_fft <= 47.55238723754883 ){
                                    return 9;
                                }else{
                                    return 0;
                                }
                            }
                        }else{
                            if ( std_ <= 0.18374106287956238 ){
                                return 9;
                            }else{
                                if ( dist_min_mean_dados <= 1.23941171169281 ){
                                    return 9;
                                }else{
                                    if ( max_aux <= 84.0 ){
                                        return 4;
                                    }else{
                                        return 5;
                                    }
                                }
                            }
                        }
                    }else{
                        if ( dist_max_mean_dados <= 1.3581883907318115 ){
                            if ( dist_max_mean_dados <= 1.0673907995224 ){
                                return 9;
                            }else{
                                if ( soma_fft <= 72.14035034179688 ){
                                    return 5;
                                }else{
                                    if ( dist_min_mean_dados <= 1.179282307624817 ){
                                        if ( max_env_dados <= 1.176605463027954 ){
                                            if ( dist_max_mean_dados <= 1.1412479877471924 ){
                                                return 4;
                                            }else{
                                                return 8;
                                            }
                                        }else{
                                            if ( min_env_dados <= 0.9835638403892517 ){
                                                return 9;
                                            }else{
                                                return 5;
                                            }
                                        }
                                    }else{
                                        if ( std_ <= 0.2146868109703064 ){
                                            return 5;
                                        }else{
                                            if ( soma_fft <= 75.38710021972656 ){
                                                return 9;
                                            }else{
                                                return 5;
                                            }
                                        }
                                    }
                                }
                            }
                        }else{
                            if ( dist_min_mean_dados <= 1.2666449546813965 ){
                                if ( min_env_dados <= 0.9804352521896362 ){
                                    return 9;
                                }else{
                                    return 4;
                                }
                            }else{
                                return 5;
                            }
                        }
                    }
                }else{
                    if ( std_ <= 0.2199658453464508 ){
                        if ( min_env_dados <= 0.9404390454292297 ){
                            return 9;
                        }else{
                            if ( soma_fft <= 35.422080993652344 ){
                                return 5;
                            }else{
                                if ( dist_min_mean_dados <= 1.2987964153289795 ){
                                    return 6;
                                }else{
                                    return 4;
                                }
                            }
                        }
                    }else{
                        if ( rms_dados <= 0.6961765289306641 ){
                            return 4;
                        }else{
                            return 9;
                        }
                    }
                }
            }
        }
    }else{
        if ( dist_env_fft_mx <= 1.868314266204834 ){
            if ( rms_dados <= 0.724662184715271 ){
                if ( dist_env_fft_mx <= 0.8871480226516724 ){
                    if ( std_ <= 0.21028323471546173 ){
                        if ( std_ <= 0.19942015409469604 ){
                            if ( dist_env_fft_mx <= 0.5835474133491516 ){
                                if ( std_ <= 0.19462385773658752 ){
                                    if ( dist_max_mean_dados <= 1.162750244140625 ){
                                        if ( std_ <= 0.19231544435024261 ){
                                            if ( max_env_dados <= 1.1520411968231201 ){
                                                return 4;
                                            }else{
                                                return 0;
                                            }
                                        }else{
                                            if ( dist_cont <= -0.5 ){
                                                return 9;
                                            }else{
                                                return 2;
                                            }
                                        }
                                    }else{
                                        if ( max_ <= 625662.25 ){
                                            return 4;
                                        }else{
                                            if ( min_env_dados <= 1.088700532913208 ){
                                                return 9;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }
                                }else{
                                    if ( min_env_dados <= 1.0315974950790405 ){
                                        if ( dist_min_mean_dados <= 1.2008538246154785 ){
                                            if ( max_ <= 177422.03125 ){
                                                return 8;
                                            }else{
                                                return 9;
                                            }
                                        }else{
                                            return 5;
                                        }
                                    }else{
                                        if ( min_env_dados <= 1.068652629852295 ){
                                            if ( max_aux <= 156.5 ){
                                                return 9;
                                            }else{
                                                return 7;
                                            }
                                        }else{
                                            if ( soma_fft <= 75.50495910644531 ){
                                                return 8;
                                            }else{
                                                return 8;
                                            }
                                        }
                                    }
                                }
                            }else{
                                if ( soma_fft <= 40.527259826660156 ){
                                    if ( rms_dados <= 0.707614541053772 ){
                                        if ( max_ <= 0.537262499332428 ){
                                            if ( rms_dados <= 0.6979764103889465 ){
                                                return 0;
                                            }else{
                                                return 9;
                                            }
                                        }else{
                                            if ( rms_dados <= 0.6891309022903442 ){
                                                return 4;
                                            }else{
                                                return 8;
                                            }
                                        }
                                    }else{
                                        if ( dist_env_fft_mx <= 0.7392328381538391 ){
                                            return 4;
                                        }else{
                                            return 3;
                                        }
                                    }
                                }else{
                                    if ( dist_max_mean_dados <= 1.1705820560455322 ){
                                        if ( soma_fft <= 82.94120788574219 ){
                                            if ( max_ <= 1915157.0 ){
                                                return 0;
                                            }else{
                                                return 8;
                                            }
                                        }else{
                                            if ( dist_env_fft_mx <= 0.6363586187362671 ){
                                                return 8;
                                            }else{
                                                return 0;
                                            }
                                        }
                                    }else{
                                        if ( soma_fft <= 74.85582733154297 ){
                                            if ( std_ <= 0.1946454644203186 ){
                                                return 9;
                                            }else{
                                                return 7;
                                            }
                                        }else{
                                            if ( dist_min_mean_dados <= 1.2705918550491333 ){
                                                return 8;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }
                                }
                            }
                        }else{
                            if ( dist_min_mean_dados <= 1.0087964534759521 ){
                                if ( max_env_dados <= 1.0207200050354004 ){
                                    return 0;
                                }else{
                                    return 7;
                                }
                            }else{
                                if ( dist_max_mean_dados <= 1.2145659923553467 ){
                                    if ( max_ <= 76595.578125 ){
                                        if ( std_ <= 0.20282770693302155 ){
                                            if ( max_aux <= 8.5 ){
                                                return 7;
                                            }else{
                                                return 7;
                                            }
                                        }else{
                                            if ( soma_fft <= 29.86145782470703 ){
                                                return 0;
                                            }else{
                                                return 8;
                                            }
                                        }
                                    }else{
                                        if ( dist_env_fft_mx <= 0.5896306037902832 ){
                                            if ( soma_fft <= 30.815536499023438 ){
                                                return 0;
                                            }else{
                                                return 7;
                                            }
                                        }else{
                                            if ( dist_max_mean_dados <= 1.0491598844528198 ){
                                                return 8;
                                            }else{
                                                return 8;
                                            }
                                        }
                                    }
                                }else{
                                    if ( dist_min_mean_dados <= 1.2124030590057373 ){
                                        if ( dist_env_fft_mx <= 0.6853326559066772 ){
                                            if ( max_env_dados <= 1.2513031959533691 ){
                                                return 7;
                                            }else{
                                                return 8;
                                            }
                                        }else{
                                            if ( soma_fft <= 83.79144287109375 ){
                                                return 8;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }else{
                                        if ( min_env_dados <= 1.0385658740997314 ){
                                            if ( dist_env_fft_mx <= 0.7214072942733765 ){
                                                return 5;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            if ( dist_env_fft_mx <= 0.581350564956665 ){
                                                return 7;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }else{
                        if ( soma_fft <= 29.3249568939209 ){
                            return 0;
                        }else{
                            if ( rms_dados <= 0.7123234272003174 ){
                                if ( dist_env_fft_mx <= 0.5920982956886292 ){
                                    if ( std_ <= 0.22103148698806763 ){
                                        if ( rms_dados <= 0.705919086933136 ){
                                            if ( rms_dados <= 0.6952102184295654 ){
                                                return 2;
                                            }else{
                                                return 8;
                                            }
                                        }else{
                                            if ( soma_fft <= 71.29129791259766 ){
                                                return 0;
                                            }else{
                                                return 8;
                                            }
                                        }
                                    }else{
                                        if ( dist_min_mean_dados <= 1.1030972003936768 ){
                                            if ( max_aux <= 132.5 ){
                                                return 9;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            return 9;
                                        }
                                    }
                                }else{
                                    if ( std_ <= 0.22460174560546875 ){
                                        if ( max_aux <= 150.5 ){
                                            if ( soma_fft <= 33.01667404174805 ){
                                                return 9;
                                            }else{
                                                return 8;
                                            }
                                        }else{
                                            if ( dist_max_mean_dados <= 1.1337635517120361 ){
                                                return 8;
                                            }else{
                                                return 8;
                                            }
                                        }
                                    }else{
                                        if ( dist_env_fft_mx_mean <= 0.5091851949691772 ){
                                            return 9;
                                        }else{
                                            if ( max_env_dados <= 1.0833890438079834 ){
                                                return 4;
                                            }else{
                                                return 0;
                                            }
                                        }
                                    }
                                }
                            }else{
                                if ( max_env_dados <= 1.2181668281555176 ){
                                    if ( std_ <= 0.2173984944820404 ){
                                        if ( dist_env_fft_mx_mean <= 0.3798717260360718 ){
                                            if ( std_ <= 0.21203967928886414 ){
                                                return 8;
                                            }else{
                                                return 0;
                                            }
                                        }else{
                                            if ( max_aux <= 133.0 ){
                                                return 3;
                                            }else{
                                                return 8;
                                            }
                                        }
                                    }else{
                                        if ( soma_fft <= 35.67333984375 ){
                                            if ( dist_max_mean_dados <= 1.0822291374206543 ){
                                                return 0;
                                            }else{
                                                return 3;
                                            }
                                        }else{
                                            if ( std_ <= 0.2394019514322281 ){
                                                return 0;
                                            }else{
                                                return 0;
                                            }
                                        }
                                    }
                                }else{
                                    if ( dist_env_fft_mx <= 0.481522798538208 ){
                                        if ( dist_max_mean_dados <= 1.225484848022461 ){
                                            return 9;
                                        }else{
                                            if ( std_ <= 0.21337756514549255 ){
                                                return 5;
                                            }else{
                                                return 3;
                                            }
                                        }
                                    }else{
                                        if ( max_env_dados <= 1.2632615566253662 ){
                                            if ( soma_fft <= 76.130126953125 ){
                                                return 8;
                                            }else{
                                                return 0;
                                            }
                                        }else{
                                            if ( std_ <= 0.2160070538520813 ){
                                                return 8;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }else{
                    if ( std_ <= 0.1785264015197754 ){
                        if ( dist_max_mean_dados <= 1.1945525407791138 ){
                            if ( rms_dados <= 0.6958750486373901 ){
                                if ( soma_fft <= 39.61693572998047 ){
                                    return 4;
                                }else{
                                    return 2;
                                }
                            }else{
                                if ( soma_fft <= 31.974624633789062 ){
                                    return 9;
                                }else{
                                    return 0;
                                }
                            }
                        }else{
                            if ( dist_max_mean_dados <= 1.292365312576294 ){
                                if ( max_aux <= 77.5 ){
                                    return 4;
                                }else{
                                    return 0;
                                }
                            }else{
                                return 9;
                            }
                        }
                    }else{
                        if ( std_ <= 0.22264182567596436 ){
                            if ( soma_fft <= 47.368656158447266 ){
                                if ( dist_min_mean_dados <= 1.0968074798583984 ){
                                    if ( max_ <= 0.01143163163214922 ){
                                        return 0;
                                    }else{
                                        if ( dist_min_mean_dados <= 1.010157823562622 ){
                                            if ( dist_env_fft_mx <= 0.9350285530090332 ){
                                                return 8;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            if ( std_ <= 0.19898205995559692 ){
                                                return 8;
                                            }else{
                                                return 8;
                                            }
                                        }
                                    }
                                }else{
                                    if ( rms_dados <= 0.7102698087692261 ){
                                        return 4;
                                    }else{
                                        if ( dist_max_mean_dados <= 1.1899938583374023 ){
                                            if ( max_ <= 478469.9375 ){
                                                return 3;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            return 4;
                                        }
                                    }
                                }
                            }else{
                                if ( std_ <= 0.1975899189710617 ){
                                    if ( std_ <= 0.1931760609149933 ){
                                        if ( rms_dados <= 0.6960960626602173 ){
                                            if ( dist_max_mean_dados <= 1.1403409242630005 ){
                                                return 8;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            if ( max_env_dados <= 1.1970874071121216 ){
                                                return 0;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }else{
                                        if ( min_env_dados <= 1.0046498775482178 ){
                                            return 4;
                                        }else{
                                            if ( max_aux <= 162.5 ){
                                                return 8;
                                            }else{
                                                return 2;
                                            }
                                        }
                                    }
                                }else{
                                    if ( std_ <= 0.21225523948669434 ){
                                        if ( rms_dados <= 0.6987156867980957 ){
                                            if ( dist_max_mean_dados <= 1.1513729095458984 ){
                                                return 8;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            if ( dist_env_fft_mx <= 1.0395660400390625 ){
                                                return 4;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }else{
                                        if ( dist_min_mean_dados <= 1.2054483890533447 ){
                                            if ( dist_env_fft_mx_mean <= 0.7514898777008057 ){
                                                return 0;
                                            }else{
                                                return 8;
                                            }
                                        }else{
                                            if ( dist_max_mean_dados <= 1.2729393243789673 ){
                                                return 3;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }
                                }
                            }
                        }else{
                            if ( rms_dados <= 0.7071486115455627 ){
                                if ( std_ <= 0.22460800409317017 ){
                                    return 8;
                                }else{
                                    if ( dist_max_mean_dados <= 1.1449341773986816 ){
                                        return 4;
                                    }else{
                                        return 8;
                                    }
                                }
                            }else{
                                if ( dist_max_mean_dados <= 1.1712682247161865 ){
                                    if ( dist_env_fft_mx_mean <= 1.4838929176330566 ){
                                        return 0;
                                    }else{
                                        return 8;
                                    }
                                }else{
                                    if ( max_env_dados <= 1.214487075805664 ){
                                        return 3;
                                    }else{
                                        return 4;
                                    }
                                }
                            }
                        }
                    }
                }
            }else{
                if ( dist_max_mean_dados <= 1.344738245010376 ){
                    if ( min_env_dados <= 1.0040428638458252 ){
                        if ( dist_env_fft_mx <= 0.8492268323898315 ){
                            if ( dist_env_fft_mx_mean <= 0.48983097076416016 ){
                                if ( std_ <= 0.22613003849983215 ){
                                    return 5;
                                }else{
                                    return 9;
                                }
                            }else{
                                if ( max_ <= 3623.5712890625 ){
                                    if ( dist_max_mean_dados <= 1.1737217903137207 ){
                                        if ( std_ <= 0.2058463990688324 ){
                                            return 3;
                                        }else{
                                            return 5;
                                        }
                                    }else{
                                        if ( dist_max_mean_dados <= 1.300421953201294 ){
                                            return 3;
                                        }else{
                                            return 4;
                                        }
                                    }
                                }else{
                                    return 3;
                                }
                            }
                        }else{
                            if ( dist_min_mean_dados <= 1.0835309028625488 ){
                                if ( soma_fft <= 31.79949378967285 ){
                                    return 3;
                                }else{
                                    if ( dist_env_fft_mx <= 1.0155739784240723 ){
                                        if ( max_env_dados <= 1.3141486644744873 ){
                                            return 4;
                                        }else{
                                            return 3;
                                        }
                                    }else{
                                        return 4;
                                    }
                                }
                            }else{
                                if ( rms_dados <= 0.7450553178787231 ){
                                    if ( dist_max_mean_dados <= 1.2243332862854004 ){
                                        return 3;
                                    }else{
                                        if ( max_ <= 5972.3115234375 ){
                                            return 3;
                                        }else{
                                            if ( min_env_dados <= 1.0003604888916016 ){
                                                return 3;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }
                                }else{
                                    return 3;
                                }
                            }
                        }
                    }else{
                        if ( max_ <= 12093.962890625 ){
                            if ( rms_dados <= 0.8056789636611938 ){
                                return 5;
                            }else{
                                return 3;
                            }
                        }else{
                            if ( dist_env_fft_mx <= 0.5945479869842529 ){
                                if ( dist_env_fft_mx <= 0.5045202374458313 ){
                                    if ( dist_max_mean_dados <= 1.171065092086792 ){
                                        if ( dist_env_fft_mx <= 0.35652682185173035 ){
                                            if ( max_env_dados <= 1.2131319046020508 ){
                                                return 3;
                                            }else{
                                                return 9;
                                            }
                                        }else{
                                            return 5;
                                        }
                                    }else{
                                        if ( max_aux <= 167.5 ){
                                            if ( dist_env_fft_mx <= 0.26699522137641907 ){
                                                return 5;
                                            }else{
                                                return 5;
                                            }
                                        }else{
                                            return 3;
                                        }
                                    }
                                }else{
                                    if ( min_env_dados <= 1.0264334678649902 ){
                                        if ( max_ <= 1169354.625 ){
                                            return 3;
                                        }else{
                                            return 5;
                                        }
                                    }else{
                                        if ( dist_min_mean_dados <= 1.2114546298980713 ){
                                            if ( rms_dados <= 0.7476444244384766 ){
                                                return 3;
                                            }else{
                                                return 5;
                                            }
                                        }else{
                                            if ( max_env_dados <= 1.3431631326675415 ){
                                                return 3;
                                            }else{
                                                return 4;
                                            }
                                        }
                                    }
                                }
                            }else{
                                if ( rms_dados <= 0.7322478294372559 ){
                                    if ( dist_max_mean_dados <= 1.188057780265808 ){
                                        return 0;
                                    }else{
                                        if ( dist_env_fft_mx <= 0.8632107973098755 ){
                                            return 3;
                                        }else{
                                            return 4;
                                        }
                                    }
                                }else{
                                    if ( min_env_dados <= 1.0608654022216797 ){
                                        if ( std_ <= 0.2076788991689682 ){
                                            if ( dist_max_mean_dados <= 1.3009002208709717 ){
                                                return 3;
                                            }else{
                                                return 4;
                                            }
                                        }else{
                                            return 3;
                                        }
                                    }else{
                                        if ( max_aux <= 146.0 ){
                                            return 5;
                                        }else{
                                            if ( max_aux <= 155.5 ){
                                                return 3;
                                            }else{
                                                return 3;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }else{
                    if ( rms_dados <= 0.7397271394729614 ){
                        if ( max_env_dados <= 1.4133496284484863 ){
                            if ( dist_env_fft_mx <= 0.49805283546447754 ){
                                return 3;
                            }else{
                                return 4;
                            }
                        }else{
                            if ( dist_env_fft_mx <= 0.6222473382949829 ){
                                if ( rms_dados <= 0.7287666201591492 ){
                                    return 5;
                                }else{
                                    return 4;
                                }
                            }else{
                                if ( std_ <= 0.2107650637626648 ){
                                    if ( max_ <= 1347416.25 ){
                                        return 8;
                                    }else{
                                        return 4;
                                    }
                                }else{
                                    return 8;
                                }
                            }
                        }
                    }else{
                        if ( dist_env_fft_mx_mean <= 0.3339883089065552 ){
                            if ( dist_max_mean_dados <= 1.390744686126709 ){
                                if ( rms_dados <= 0.762340784072876 ){
                                    return 3;
                                }else{
                                    if ( rms_dados <= 0.807884693145752 ){
                                        if ( max_aux <= 163.0 ){
                                            if ( max_env_dados <= 1.3969590663909912 ){
                                                return 5;
                                            }else{
                                                return 3;
                                            }
                                        }else{
                                            return 3;
                                        }
                                    }else{
                                        return 3;
                                    }
                                }
                            }else{
                                return 3;
                            }
                        }else{
                            if ( rms_dados <= 0.7765715718269348 ){
                                if ( dist_min_mean_dados <= 1.2156901359558105 ){
                                    if ( soma_fft <= 39.669090270996094 ){
                                        return 3;
                                    }else{
                                        return 4;
                                    }
                                }else{
                                    if ( rms_dados <= 0.7511482238769531 ){
                                        if ( std_ <= 0.23109543323516846 ){
                                            if ( min_env_dados <= 1.031557321548462 ){
                                                return 4;
                                            }else{
                                                return 3;
                                            }
                                        }else{
                                            return 3;
                                        }
                                    }else{
                                        return 3;
                                    }
                                }
                            }else{
                                return 3;
                            }
                        }
                    }
                }
            }
        }else{
            if ( dist_env_fft_mx <= 4.34155797958374 ){
                if ( rms_dados <= 0.7719789743423462 ){
                    if ( min_env_dados <= 1.0377964973449707 ){
                        if ( dist_max_mean_dados <= 1.1519994735717773 ){
                            if ( soma_fft <= 56.571197509765625 ){
                                if ( dist_max_mean_dados <= 1.0239150524139404 ){
                                    if ( rms_dados <= 0.7007275819778442 ){
                                        return 4;
                                    }else{
                                        if ( dist_env_fft_mx <= 2.1926703453063965 ){
                                            return 4;
                                        }else{
                                            return 6;
                                        }
                                    }
                                }else{
                                    return 4;
                                }
                            }else{
                                if ( max_aux <= 151.0 ){
                                    if ( std_ <= 0.18255814909934998 ){
                                        return 2;
                                    }else{
                                        return 6;
                                    }
                                }else{
                                    if ( soma_fft <= 86.91825866699219 ){
                                        return 4;
                                    }else{
                                        return 1;
                                    }
                                }
                            }
                        }else{
                            if ( dist_env_fft_mx <= 1.8774800300598145 ){
                                if ( max_aux <= 5.0 ){
                                    return 8;
                                }else{
                                    return 4;
                                }
                            }else{
                                return 4;
                            }
                        }
                    }else{
                        if ( max_aux <= 138.5 ){
                            return 6;
                        }else{
                            if ( dist_max_mean_dados <= 1.2679085731506348 ){
                                if ( soma_fft <= 84.1214599609375 ){
                                    if ( std_ <= 0.1963728368282318 ){
                                        return 4;
                                    }else{
                                        return 6;
                                    }
                                }else{
                                    if ( min_env_dados <= 1.1217641830444336 ){
                                        if ( dist_env_fft_mx <= 3.175898313522339 ){
                                            return 4;
                                        }else{
                                            return 6;
                                        }
                                    }else{
                                        return 6;
                                    }
                                }
                            }else{
                                if ( std_ <= 0.21854904294013977 ){
                                    return 4;
                                }else{
                                    return 8;
                                }
                            }
                        }
                    }
                }else{
                    return 3;
                }
            }else{
                if ( min_env_dados <= 1.0048415660858154 ){
                    return 4;
                }else{
                    if ( min_env_dados <= 1.0985643863677979 ){
                        if ( max_env_dados <= 1.3858789205551147 ){
                            return 6;
                        }else{
                            return 4;
                        }
                    }else{
                        return 6;
                    }
                }
            }
        }
    }
}

void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
  sampling_period_us = round(1000000*(1.0/SAMPLING_FREQUENCY));
}

String ler(){
  String conteudo = "";
  char caractere;
  while(conteudo == ""){
    while(Serial.available() > 0) {
      // Lê byte da serial
      caractere = Serial.read();
      // Ignora caractere de quebra de linha
      if (caractere != '\n'){
        // Concatena valores
        conteudo.concat(caractere);
      }
      // Aguarda buffer serial ler próximo caractere
      delay(10);
    }
  }
  return conteudo;
}

float dados[SAMPLES_DATA];
float f[12];

void loop() {
  
  for (int i=0; i<SAMPLES_DATA; i++){
    dados[i]=ler().toFloat(); 
    Serial.println(i);
  }

  unsigned long StartTime = millis();
  construct_arg(f, dados);
  
  int classe = classifier(f[0], f[1], f[2], f[3], f[4], f[5],
        f[6], f[7], f[8], f[9], f[10], f[11]);
  
  Serial.println(classe);
  
  unsigned long CurrentTime = millis();
  unsigned long ElapsedTime = CurrentTime - StartTime;
  Serial.println(ElapsedTime);
}
