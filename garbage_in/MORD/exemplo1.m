%Exempl oda aplica��o de um detector univariado 

clear all, close all, clc

%Gera��o do sinal sint�tico ---------------------------------------
tj = 256; %tamanho de cada janela 
M = 50;  %n�mero da janela
fs = 256;  %frequ�ncia de amostragem
y = randn(tj*M,1); %ruido sint�tico 
f = 12; %frequencia do estimulo 
A = 1; %amplitude do sinal 
aux = A*cos(2*pi*f*[1:length(y)]/tj)';
y = y+A*cos(2*pi*f*[1:length(y)]/tj)';

% figure
% subplot(211);plot(aux);xlim([0 400])
% title('Sinal senoidal');ylabel('Amplitude')
% subplot(212);plot(y);xlim([0 400])
% title('Sinal EEG sint�tico');xlabel('tempo (s)');ylabel('Amplitude')
%-------------------------------------------------------------------

%detector
detector = 'MSC'; %detector univariado: 'MSC', 'CSM', 'LFT'


%parametros -------------------------------------------------
parameters.tj = tj;
parameters.M = M;
parameters.N = 1;
parameters.fs = fs;
%-----------------------------------------------------------

[value] = mord(detector,y,parameters); %obter o valor do detector 

%valor cr�tico ------------------------------------------------
alpha = 0.05; %n�vel de signific�ncia
methods = 'theoretical'; %m�todo para obter o valor cr�tico
CV  = critical_value(alpha,detector,parameters,methods) ;

%imprimir ---------------------------------------------
f = 0:parameters.fs/parameters.tj:parameters.fs/2;
figure 
plot(f,value,'k')
hold on 
plot([0 max(f)],[CV CV],'r')
xlabel('frequ�ncia','fontsize',12)
ylabel('Valor da ORD','fontsize',12) 
legend({'ORD','CV'})




