%exemplo da utilização da função MORD 
clear all, close all, clc 
% addpath  MORD  

%Definir os parametros ---------------------------------------------

parameters.M=50; %número de janelas 
parameters.tj=256; %número de pontos em cada janela 
parameters.N= 2 ; %número de canais 
parameters.fs= 256 ; %frequencia de amostragem 

%Gerar o sinal multivariado ------------------------------------------
x = randn(parameters.M*parameters.tj,parameters.N);
f = 12; 
x = x+repmat(sin(2*pi*f*[1:parameters.M*parameters.tj]/parameters.tj)',1,parameters.N);
%------------------------------------------------------------------------------


detector = 'MMSC';  %definir o detector multivariado 
%Detector = {'MSC', 'CSM', 'LFT', 'MMSC', 'aMSC', 'pMSC', 'MCSM', 'aCSM', 'pCSM', 'MLFT', 'aLFT',...



value = mord(detector,x,parameters); %gera o valor do detector 
%observação: parameters é uma estrutura com os parametros utilizados no
%%detector (string): nome do detector 
%parameters: parametros do detector 
% parameters.M -> numero de janelas 
% parameters.L -> banda Lateral
% parameters.N -> número de canais
%methods -> método utilizado para estimar o valor crítico
%   pararametrs.
%yin, xin -> sinais utilizadas para as correções (cholesky correção)

%% Geração do valor crítico --------------------------------------------

alpha = 0.05;
methods = 'time_Cholesky_corrected'; %método utilizado para gerar o valor crítico
%list_methods = {'theoretical', 'Monte_Carlo_default','time_Cholesky_corrected','frequency_Cholesky_corrected'};
parameters.mistura = 'cholesky'; %no caso do método 'time_Cholesky_corrected' foi utilizado o método de cholesky obter a matriz de mistura
CV  = critical_value(alpha,detector,parameters,methods,x) ;
%-----------------------------------------------------------------------------
%% 
f = 0:parameters.fs/parameters.tj:parameters.fs/2;
figure 
plot(f,value,'k')
hold on 
plot([0 max(f)],[CV CV],'r')
xlabel('frequencia','fontsize',12)
ylabel('Valor da MORD','fontsize',12) 
legend({'MORD','CV'})




