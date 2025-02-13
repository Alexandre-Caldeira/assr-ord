%% Traçando graficos
% close all;
clearvars;clc;
% Carrega variavies (têm que estar na pasta local do caminho)
load('Workspace_M1','FP','TXD', 'timeM');
fp1 = FP;
txd1 = TXD;
timeM1 = timeM;
load('Workspace_M2','FP','TXD', 'timeM');
fp2 = FP;
txd2 = TXD;
timeM2 = timeM;
load('Workspace_M3','FP','TXD', 'timeM');
fp3 = FP;
txd3 = TXD;
timeM3 = timeM;
%%
load('Workspace_M4','FP','TXD', 'timeM');
fp4 = FP;
txd4 = TXD;
timeM4 = timeM;
load('Workspace_M1','linf');
load('Workspace_M1','lsup');
%% Teste McNemar
% 50dB Single test first line
%Comparar Método 04 x Standart
%Comparar Método 04 x Método 02
%% Tempo
t4 = min(timeM4); tm4 = mean(timeM4);
t3 = min(timeM3); tm3 = mean(timeM3);
t2 = min(timeM2); tm2 = mean(timeM2);
t1 = min(timeM1); tm1 = mean(timeM1);
t = timeM4(1);
menort4 = (t4-t)/t;
menort3 = (t3-t)/t;
menort12 = (t2-t)/t;
%% Taxa de Detecção
dr4 = max(txd4);
dr3 = max(txd3);
dr2 = max(txd2);
dr1 = max(txd1);
dr = txd4(1);
maiordr4 = (dr4-dr)/dr;
maiordr3 = (dr3-dr)/dr;
maiordr12 = (dr2-dr)/dr;
%%
% %% Pontos Fora do Intervalo - M1
% %Pontos acima
% Ac1 = mean(fp1*100 > lsup);
% %Pontos abaixo
% Ab1 = mean(fp1*100 < linf); 
% %% Pontos Fora do Intervalo - M2
% %Pontos acima
% Ac2 = mean(fp2*100 > lsup);
% %Pontos abaixo
% Ab2 = mean(fp2*100 < linf); 
% %% Pontos Fora do Intervalo - M3
% %Pontos acima
% Ac3 = mean(fp3*100 > lsup);
% %Pontos abaixo
% Ab3 = mean(fp3*100 < linf); 
% %% Pontos Fora do Intervalo - M4
% %Pontos acima
% Ac4 = mean(fp4*100 > lsup);
% %Pontos abaixo
% Ab4 = mean(fp4*100 < linf); 

% load TXD_correcao
% load TXD_sem_correcao
% 
% aumento=TXD_correcao-TXD_sem_correcao;
% for ii=1:size(TXD_correcao,1)   
%         aumenta_TXD(ii)=aumento(ii)/TXD_sem_correcao(ii);  
% end   
% aumenta_TXD=aumenta_TXD*100;
% 
% mean(aumenta_TXD)
% 
% median(aumenta_TXD)
% 
% x=[1:1:1342];
% figure
% plot(x,aumenta_TXD,'Linewidth',1.2)
% xlabel('Indice do Protocolo','Fontsize',12)
% ylabel('Aumento da Taxa de detecção(%)','Fontsize',12)
