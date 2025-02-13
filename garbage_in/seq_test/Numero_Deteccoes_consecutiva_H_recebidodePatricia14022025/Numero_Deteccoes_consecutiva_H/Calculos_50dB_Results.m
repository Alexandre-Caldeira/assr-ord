%% Traçando graficos 30dB
% close all;
clearvars;clc;
%% Carrega variavies (têm que estar na pasta local do caminho)
%  01
load('Workspace_M1','FP','TXD', 'timeM', 'Tdr', 'Ttime');
fp1 = FP; txd1 = TXD; timeM1 = timeM; 
Tdr1 = Tdr;
Ttime1 = Ttime;
%% 02
load('Workspace_M2','FP','TXD', 'timeM', 'Tdr');
fp2 = FP;
txd2 = TXD;
timeM2 = timeM;
Tdr2 = Tdr;
%% 03
load('Workspace_M3','FP','TXD', 'timeM', 'Tdr');
fp3 = FP;
txd3 = TXD;
timeM3 = timeM;
Tdr3 = Tdr;
%% 04
load('Workspace_M4','FP','TXD', 'timeM', 'Tdr');
fp4 = FP;
txd4 = TXD;
timeM4 = timeM;
Tdr4 = Tdr;
load('Workspace_M1','linf', 'binsM');
load('Workspace_M1','lsup');
%% 05
load('Workspace_NND2_M4','FP','TXD', 'timeM', 'Tdr');
fp4 = FP;
txd4 = TXD;
timeM4 = timeM;
Tdr4 = Tdr;
load('Workspace_M1','linf', 'binsM');
load('Workspace_M1','lsup');



%% Teste McNemar
% h = testcholdout(YHat1,YHat2,Y) returns the test decision, by conducting the mid-p-value McNemar test, 
% from testing the null hypothesis that the predicted class labels YHat1 and YHat2 have equal accuracy for
%  predicting the true class labels Y. The alternative hypothesis is that the labels have unequal accuracy.
% h = 1 indicates to reject the null hypothesis at the 5% significance level. 
% h = 0 indicates to not reject the null hypothesis at 5% level.
% 30dB Single test last line
%Comparar Método 04 x Standart
%Comparar Método 04 x Método 02
% freqEstim = [81 83 85 87 89 91 93 95]; %frequencia de estimulação fixa 
% binsM = [freqEstim] +1;%bins de modulação
%Tdr = 120 x 1342 x 11
% [valores, indices] = min(Tdr4); % devolve valor e indices em Tdr dos menores tempos de detecção
% Ppareto = P(indices);

% Tdr2 = reshape( Tdr1, [88,1]);
%% Descobrir parâmetros de Pareto
%Indices paretos: function: [p, idxs] = paretoFront([TXD,(-timeM)] );
%Método 1 [489;409;1268;1127;503;545;646;1]     - 1
%Método 2 [489;409;1268;1127;598;799;646;362;8] - 362
%Método 3 [1092;484;1074;808;566;113;174;439;1] - 1
%Método 4 [1338;1223;1083; 1092;484;483;1074;808;566; 113; 174;1] - 1
[p, idxs] = paretoFront([TXD,(-timeM)] );
%%
p1 =  1; p2 =  362; p3 =  1; p4 =  1;  % 67,5% 
% p1 =  1; p2 =  8; p3 =  1; p4 =  1;  % últimaª 
% p1 =  489; p2 =  489; p3 =  1092; p4 =  1338; % 1ª 
% p41 = 1338; p42 = 1223; p43 = 1083; % 3 1ª M4 
p43 =  113; p44 = 1; % 3ª e 4ª M4 
% p43 = 174; p4 = 1;  %11ª e 12ª M4
% p43 = 174; p4 = 1;  %4ª e 11ª M4 
%%
Vdr4 = Tdr4 ( binsM, p44 ,:);  % Selecionando freq Est. Melhor P voluntários
Vdr4 = Vdr4(:); 
Vdr3 = Tdr4 ( binsM, p43 ,:);  % Selecionando freq Est. Melhor P voluntários
Vdr3 = Vdr3(:); 
h4 = testcholdout (Vdr4, Vdr3, V);
%%
Vdr4 = Tdr4 ( binsM, p44 ,:);  % Selecionando freq Est. Melhor P voluntários
Vdr4 = Vdr4(:); 
Vdr3 = Tdr4 ( binsM, 1074 ,:);  % Selecionando freq Est. Melhor P voluntários
Vdr3 = Vdr3(:);
Vdrs = Tdr4 (binsM, 1, :); % Vetor single test last line
Vdrs = Vdrs(:);             
h4 = testcholdout (Vdr4, Vdr3, V);
%%
Vdr4 = Tdr4 ( binsM, p4 ,:);  % Selecionando freq Est. Melhor P voluntários
Vdr4 = Vdr4(:);             % Vetorizando
Vdr3 = Tdr3 ( binsM, p3 ,:);  % Selecionando freq Est. Melhor P voluntários
Vdr3 = Vdr3(:);             % Vetorizando
Vdr2 = Tdr2 ( binsM, p2 ,:);  % Selecionando freq Est. Melhor P voluntários
Vdr2 = Vdr2(:);             % Vetorizando
Vdr1 = Tdr1 ( binsM, p1 ,:);  % Selecionando freq Est. Melhor P voluntários
Vdr1 = Vdr1(:);             % Vetorizando
Vdrs = Tdr1 (binsM, 1, :); % Vetor single test last line
Vdrs = Vdrs(:);             % Vetorizando
V = ones (size(Vdr1));
% h = testcholdout (Vdr1, Vdrs, V);

h1 = testcholdout (Vdr1, Vdrs, V);
h2 = testcholdout (Vdr2, Vdrs, V);
h3 = testcholdout (Vdr3, Vdrs, V);
h4 = testcholdout (Vdr4, Vdrs, V);

% h = testcholdout (Vdr4, Vdr1, V);

%% Tempo
t4 = min(timeM4); tm4 = mean(timeM4);
t3 = min(timeM3); tm3 = mean(timeM3);
t2 = min(timeM2); tm2 = mean(timeM2);
t1 = min(timeM1); tm1 = mean(timeM1);
t = timeM4(1); %single test first line
mediat4 = (tm4-t)/t;
menort4 = (t4-t)/t;
menort3 = (t3-t)/t;
menort12 = (t2-t)/t; %tempo min 1 e 2 são =
%% Taxa de Detecção
dr4 = max(txd4);drm4 = mean(txd4);
dr3 = max(txd3);drm3 = mean(txd3);
dr2 = max(txd2);drm2 = mean(txd2);
dr1 = max(txd1);drm1 = mean(txd1);
dr = txd4(end); %single test last line
maiordr4 = (dr4-dr)/dr;
% maiordr3 = (dr3-dr)/dr;
% maiordr2 = (dr2-dr)/dr;
maiordr1 = (dr1-dr)/dr;
