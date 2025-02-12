%fun��o com os parametros globais para todas ass simula��es 

Nruns = 10000;
parametro = 0;
alfa=.05; %significance level
M=12;% number of epochs
n_sinais = [1 2 4 6 8];% number of signals
npontos = 32; %numero de pontos de cada janela para MC
L = npontos*M; %tamanho do sinal em pontos para MC


fe = 6; 
tj = 512; %tamanho da janela para analise do EEG
fs = 256; %frequencia de amostragem do EEG
vet_arq = {'1_F6' '2_F6' '3_F6' '4_F6' '5_F6' '6_F6' '7_F6' '8_F6' '9_F6' '10_F6'};
vet_arqR = {'1_F6' '2_F6' '3_F6' '4_F6' '5_F6' '6_F6' '7_F6' '8_F6' '9_F6' '10_F6'};
vet_canais = [2,1,4,3,6,5,11,10];
nh = fix(fs/2/fe)-1;% n�mero de harmonico
binS = 1+(1:nh)*fe*tj/fs;
binR = 7+(1:nh)*fe*tj/fs;
fo = 6; 
bandL = 12; 

vetor_teste_F = [-6:6]; 
vetor_teste_R = [-9:-8,-4:-1,0,1:4,8:9];
