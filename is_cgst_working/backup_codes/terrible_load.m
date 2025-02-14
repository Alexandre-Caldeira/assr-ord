function SIGNALS = terrible_load(intens,volunt)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

cont_int = intens;
cont_vol = volunt;

% SET THIS PATH:
path = 'C:\PPGEE\SBEB_CBA_24\ASSR - Coleta OFFLINE';
addpath(path)
%vetor dos voluntários 
Vvoluntario = {'Abdon';'Ana';'BBB';'Colatina';'Erick';'Luciana';...
    'Sombra';'Quenaz';'Vinicius';'Sacola';'Wreikson'}; 

%vetor da intensidade 
Vintensidade = {'70';'60';'50';'40';'30'}; 
load('eletrodos.mat')
ganho  = 200;
remoc = 0.1/ganho; 

voluntario = cell2mat(Vvoluntario(cont_vol,:));

if cont_int==6
     load([voluntario 'ESP'], 'x','Fs','binsM','freqEstim') 
else
    intensidade = cell2mat(Vintensidade(cont_int,:));
    load([voluntario '_'  intensidade 'dB'], 'x','Fs','binsM','freqEstim') 
end

nfft = Fs;%1segundo de sinal 
 
%retirar componente DC por janela (fiz isso pq no processamento em
%tempo real é por janela)
x = x - repmat(mean(x),nfft,1); %tirar a média de cada de cada trecho - devido a remoção

% excluir os dois primeiros segundos do inicio da coleta 
x(:,1:2,:) =[]; 

% %encontrar o valor máximo por canal 
Vmax = squeeze(max(max(abs(x)),[],3));
ind = Vmax>remoc;
% xmedia = squeeze(mean(x(:,~ind,:),2));

pos_eletrodo= 1;
xmedia = x(:,~ind,pos_eletrodo);

SIGNALS = fft(xmedia,Fs);%*2/nfft*1e9;
FS = Fs;
NFFT = Fs;
SIGNALS = SIGNALS(1:floor(end/2)+1,:); % only half the FFT spectrum is valid

end