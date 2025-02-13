% COmparar detecções dos métodos
% Load das detecões de cada método
% 50db 240 janelas: 120x1342x11
clearvars
close all
clc
%% Método 01
load(['Tdr_50db_M1' '.mat'],'Tdr');
Tdr1 = Tdr (:, :, :);
%%  Método 02
load(['Tdr_50db_M2' '.mat'],'Tdr');
Tdr2 = Tdr (:, :, :);
%%  Método 03
load(['Tdr_50db_M3' '.mat'],'Tdr');
Tdr3 = Tdr (:, :, :);
%%  Método 04
% load(['Tdr_50db_M3' '.mat'],'Tdr');
% Tdr4 = Tdr (:, :, :);
%% Comparações 1 e 2
% Dois métodos detectaram -  x1*x2
Det112 = (Tdr1).*(Tdr2);
% Dois métodos não detectaram - (~x1)*(~x2)
nDet112 = (~Tdr1).*(~Tdr2);
% Método 1 sim 2 não - (x1)*(~x2)
Det12 = (Tdr1).*(~Tdr2);
% Método 2 sim 1 não - (~x1)*(x2)
Det21 = (~Tdr1).*(Tdr2);
% Tdet1 = mean (Det1');
Tdet12 = sum (sum(sum(Det12)));
Tdet21 = sum (sum(sum(Det21)));
%% Comparações 2 e 3
% Dois métodos detectaram -  x1*x2
Det123 = (Tdr2).*(Tdr3);
% Dois métodos não detectaram - (~x1)*(~x2)
nDet123 = (~Tdr2).*(~Tdr3);
% Método 1 sim 2 não - (x1)*(~x2)
Det23 = (Tdr2).*(~Tdr3);
% Método 2 sim 1 não - (~x1)*(x2)
Det32 = (~Tdr2).*(Tdr3);
% Tdet1 = mean (Det1');
Tdet2 = sum (sum(sum(Det123)));
Tdet3 = sum (sum(sum(nDet123)));
Tdet23 = sum (sum(sum(Det23)));
Tdet32 = sum (sum(sum(Det32)));
