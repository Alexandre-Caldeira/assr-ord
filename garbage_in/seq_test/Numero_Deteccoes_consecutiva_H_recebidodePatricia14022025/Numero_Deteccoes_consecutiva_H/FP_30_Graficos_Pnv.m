%% Traçando graficos
close all;clearvars;clc;
% Carrega variavies (têm que estar na pasta local do caminho)
load('Workspace_M1_30','FP');
fp1 = FP;
load('Workspace_M2_30','FP');
fp2 = FP;
load('Workspace_M3_30','FP');
fp3 = FP;
load('Workspace_M4_30','FP');
fp4 = FP;
load('Workspace_M1','linf');
load('Workspace_M1','lsup');
%%
% cria legendas
leg0 = 'FP Desejado = 5%';
legc1 = ['Confiança = ',num2str(linf),'%'];
legc2 = ['Confiança = ',num2str(lsup),'%'];
legC = ['Confiança = [',num2str(linf),' - ',num2str(lsup),'] %'];
leg1 = ['Método 1, média = ',num2str(round(mean(fp1)*100,2)),'%'];
leg2 = ['Método 2, média = ',num2str(round(mean(fp2)*100,2)),'%'];
leg3 = ['Método 3, média = ',num2str(round(mean(fp3)*100,2)),'%'];
leg4 = ['Método 4, média = ',num2str(round(mean(fp4)*100,2)),'%'];
tamanho = 10;

%% plota 1
figure

%Traços média FP e  de Limites de confiança
x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym1 = mean(fp1)*100*ones(1,length(x)); %Média Método 01
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

%Plotar Traços
plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on
plot(x,ym1,'Color','k','LineStyle',':','LineWidth',2); %Média Método 01
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

%Plotar pontos/conjuntos de FP
%Método 01
scatter(1:size(fp1,1),fp1*100,tamanho,'filled','CData',[0, 0.4470, 0.7410]);

%Legendas, títulos e Eixos
legend({leg0,leg1,legC},'fontsize',12)                %Método 01
ylabel('Falso Positivo','fontsize',14);
xlabel('Indíce dos conjuntos de Parâmetros','fontsize',14);
% title(['Intensidade de ', (Intensidade{1}), ' SPL: Métodos 01'])

%Fixar eixos
xlim([0,length(FP)])
ylim([0  9])
%% 2
figure

%Traços média FP e  de Limites de confiança
x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym2 = mean(fp2)*100*ones(1,length(x)); %Média Método 02
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

%Plotar Traços
plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on
plot(x,ym2,'Color','k','LineStyle',':','LineWidth',2);  %Média Método 02
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

%Plotar pontos/conjuntos de FP
%Método 02
scatter(1:size(fp2,1),fp2*100,tamanho,'filled','CData',[0.8500, 0.3250, 0.0980]);

%Legendas, títulos e Eixos
legend({leg0,leg2,legC},'fontsize',12)                   %Método 02
ylabel('Falso Positivo','fontsize',14);
xlabel('Indíce dos conjuntos de Parâmetros','fontsize',14);
% title ('Intensidade de 50 dB SPL');

%Fixar eixos
xlim([0,length(FP)])
ylim([0  9])
%% 3
figure

%Traços média FP e  de Limites de confiança
x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym3 = mean(fp3)*100*ones(1,length(x)); %Média Método 03
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

%Plotar Traços
plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on
plot(x,ym3,'Color','k','LineStyle',':','LineWidth',2); %Média Método 03
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

%Plotar pontos/conjuntos de FP
%Método 03
scatter(1:size(fp3,1),fp3*100,tamanho,'filled','CData',[0.9290, 0.6940, 0.1250]);

%Legendas, títulos e Eixos
legend({leg0,leg3,legC},'fontsize',12)                    %Método 03
ylabel('Falso Positivo','fontsize',14);
xlabel('Indíce dos conjuntos de Parâmetros','fontsize',14);
% title ('Intensidade de 50 dB SPL');

%Fixar eixos
xlim([0,length(FP)])
ylim([0  9])


%% 4
figure

%Traços média FP e  de Limites de confiança
x = [0 length(FP)];
y = 0.05*100*ones(1,length(x));
ym3 = mean(fp4)*100*ones(1,length(x)); %Média Método 04
ysup = lsup*ones(1,length(x));
yinf = linf*ones(1,length(x));

%Plotar Traços
plot(x,y,'Color','b','LineStyle','--','LineWidth',2);
hold on
plot(x,ym3,'Color','k','LineStyle',':','LineWidth',2);%Média Método 04
plot(x,yinf, 'Color','k','LineStyle','--','LineWidth',2);
plot(x,ysup, 'Color','k','LineStyle','--','LineWidth',2);

%Plotar pontos/conjuntos de FP 
%Método 04
scatter(1:size(fp4,1),fp4*100,tamanho,'filled','CData',[0.4940, 0.1840, 0.5560]);

%Legendas, títulos e Eixos
legend({leg0,leg4,legC},'fontsize',12)                   %Método 04
ylabel('Falso Positivo','fontsize',14);
xlabel('Indíce dos conjuntos de Parâmetros','fontsize',14);
% title ('Intensidade de 50 dB SPL');

%Fixar eixos
xlim([0,length(FP)])
ylim([0  9])