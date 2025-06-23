%% Traçando graficos
close all;clearvars;clc;
% Carrega variavies (têm que estar na pasta local do caminho)
fp1 = load('Workspace_M1','FP');
fp2 = load('Workspace_M2','FP');
fp3 = load('Workspace_M3','FP');
fp4 = load('Workspace_M4','FP');
linf = load('Workspace_M1','linf');
lsup = load('Workspace_M1','lsup');
%%
% cria legendas
leg0 = 'FP Desejado';
leg1 = ['Método 1, média = ',num2str(mean(fp1.FP))];
leg2 = ['Método 2, média = ',num2str(mean(fp2.FP))];
leg3 = ['Método 3, média = ',num2str(mean(fp3.FP))];
leg4 = ['Método 4, média = ',num2str(mean(fp4.FP))];


tamanho = 20;

%% plota 1
figure
yline(0.05*100, '-k','LineWidth',2);
ylabel('Falso Positivo','fontsize',12);
hold on
m1 = yline(mean(fp1.FP)*100,'--','LineWidth',2);
m1.Color = 'k';

p1 = scatter(1:size(fp1.FP,1),fp1.FP*100,tamanho,'filled');

p1.CData = [0, 0.4470, 0.7410];
% p1.MarkerEdgeAlpha = 0.2;

yline(linf, ':r','LineWidth',2);
yline(load('Workspace_M1','lsup').lsup, ':r','LineWidth',2);
legend(leg0,leg1);

%% 2
figure
yline(0.05*100, '-k','LineWidth',2);
ylabel('Falso Positivo','fontsize',12);
hold on
m1 = yline(mean(fp2.FP)*100,'--','LineWidth',2);
m1.Color = 'k';

p1 = scatter(1:size(fp2.FP,1),fp2.FP*100,tamanho,'filled');

p1.CData = [0.8500, 0.3250, 0.0980];
% p1.MarkerEdgeAlpha = 0.2;

yline(linf, ':r','LineWidth',2);
yline(lsup, ':r','LineWidth',2);
legend(leg0,leg2);

%% 3
figure
yline(0.05*100, '-k','LineWidth',2);
ylabel('Falso Positivo','fontsize',12);
hold on
m1 = yline(mean(fp3.FP)*100,'--','LineWidth',2);
m1.Color = 'k';

p1 = scatter(1:size(fp3.FP,1),fp3.FP*100,tamanho,'filled');

p1.CData = [0.9290, 0.6940, 0.1250];
% p1.MarkerEdgeAlpha = 0.2;

yline(load('Workspace_M3','linf').linf, ':r','LineWidth',2);
yline(load('Workspace_M3','lsup').lsup, ':r','LineWidth',2);
legend(leg0,leg3);

%% 4
figure
yline(0.05*100, '-k','LineWidth',2);
ylabel('Falso Positivo','fontsize',12);
hold on
m1 = yline(mean(fp4.FP)*100,'--','LineWidth',2);
m1.Color = 'k';

p1 = scatter(1:size(fp4.FP,1),fp4.FP*100,tamanho,'filled');

p1.CData = [0.4940, 0.1840, 0.5560];
% p1.MarkerEdgeAlpha = 0.2;

yline(load('Workspace_M3','linf').linf, ':r','LineWidth',2);
yline(load('Workspace_M3','lsup').lsup, ':r','LineWidth',2);
legend(leg0,leg4);
