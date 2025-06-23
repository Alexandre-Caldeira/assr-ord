%Aplicação do protocolo ao banco de dados 
clear all, close all, clc 

caminho = 'C:\Users\Crohma\Documents\PNV\Simulation\Testes_Seq\Sinais_EEG\';

%vetor dos voluntários 
Vvoluntario = {'Ab';'An';'Bb';'Er';'Lu';...
    'So';'Qu';'Vi';'Sa';'Ti';'Wr'}; %vetor dos voluntário 

% Vvoluntario = Vvoluntario([1,2]);

%Intensidade ----------------------
%intensidade = {'70dB';'60dB';'50dB';'40dB';'30dB';'ESP'}; %quais intensidade analisadas 
%sugestão
%vetor_Mmax = [50;50;240;440;440;20]; %número máximo de janela para cada intensidade
Intensidade = {'ESP'};
% Intensidade = {'ESP'};
% Intensidade = {'60dB'};

Mmax = 50; %valor máximo 

%% Parametros do protocolo de detecção. 

%parametros = [Min Mstep Mmax alfa_corrigido]
%parametros = [200 1 200 .05];
% nRuns  = 10000;
% Mmax = 240; %número máximo de janela 
% alfa_teste = 0.05;
% FP_desejado =0.05;
% [alfa_corrigido,NDC_minimo,cost_alfa, P] = funcao_NDC_alfaCorrigido_Mmax(nRuns,Mmax ,alfa_teste,FP_desejado);
alfa = 0.05;
FP_desejado = 0.05; 


% load(['NDC_AlfaCorrigido_MSC_Mmax' num2str(Mmax) '_alfa_'  num2str(alfa) '_FPdesejado' num2str(FP_desejado) '.mat'],'alfa_corrigido', ...
%     'NDC_minimo','P', 'nRuns')
load(['NDC_AlfaCorrigido_Mmax' num2str(Mmax) '_alfa_'  num2str(alfa) '_FPdesejado' num2str(FP_desejado) '.mat'],'alfa_corrigido', ...
    'NDC_minimo','P', 'nRuns')

parametros = [P, NDC_minimo,alfa_corrigido];
% parametros(:,5) = 0.05;


%% deixar o alfa corrigido em 0.05 nesse caso não há correçao
% for ii=1:size(parametros,1)
%   parametros(ii,5)=0.05; 
%     
% end

%%
load([caminho 'eletrodos.mat'])
%aux=[1:10];
pos_ele =1; 

ganho  = 200;
alpha = 0.05; 

%% --------------------------------------------------
remoc = [.1]/ganho; 


%% 
%******poder fazer por intensidade aqui -------
for ncanal=1:16 %eu que fiz
    for cont_vol = 1:size(Vvoluntario,1) %fazer por voluntário 
    
    voluntario = cell2mat(Vvoluntario(cont_vol,:)); %carregar o voluntário 
    intensidade = cell2mat(Intensidade); %intensidadde 
    load([caminho voluntario intensidade], 'x','Fs','binsM','freqEstim')   
  
     x = x(:,:,pos_ele); % todas as linhas e colunas e só o eletrodo utilizado %esse é o verdadeiro
     %%%%%%x=x(:,:,ncanal); %eu que alterei
     nfft = Fs;%1segundo de sinal 
         
     %retirar componente DC por janela (fiz isso pq no processamento em
     %tempo real é por janela)
     x = x - repmat(mean(x),nfft,1); %tirar a média de cada de cada trecho - devido a remoção
        
     %scluir os dois primeiros segundos do inicio da coleta 
     x(:,1:2,:) =[]; 
        
        
     %encontrar o valor máximo por canal 
      Vmax = max(abs(x),[],1);
      ind = Vmax>remoc;
      [sum(ind) cont_vol ];
      x = x(:,~ind); %removor o ruído amplitude 
      x = x(:,1:Mmax);%limitar o tamanho para o valor máximo. 
     
      %******** fazer por canal diferente ----for nCanal = 1:16 %
      
      [dr,time] = protocolo_deteccao(x, parametros);
      
      Tdr(:,:,cont_vol) = dr; %dr é=1 se atingiu ndc e zero caso contrário(dr é uma matriz onde
%as linhas representam as frequências e colunas representa o indice dos conjuntos de parametros)     
      Ttime(:,:,cont_vol) = time; %salva dr em função do voluntário
    end           
     %%fazendo a análise de desemprenho por canal
     %da qui para baixo eu que fiz a modificação
%TXD(:,ncanal) = mean(mean(Tdr(binsM,:,:),3),1)';

%Tdr(bimsm,:,:) representa o indice da frequencias de estimulação para todos os
%parametros e voluntários
TXD = mean(mean(Tdr(binsM,:,:),3),1)'; %calcula a media da detecção em função dos voluntarios 
%e em seguida em função dos parametros

binsR = binsM+1;
binsR = 1:100; 
binsR(binsM) = []; 
binsR(1:2) = []; 

FP(:,ncanal) = mean(mean(Tdr(binsR,:,:),3),1)';

timeM = time(binsM,:); 
timeM(timeM==-1) = Mmax;
timeMedio(:,ncanal) = mean(timeM,1)'*1; %1segundo por janela

end
%% Análise de desempenho 

%TXD - analisar as freq. estimulação 
%binsM = [82    84    86    88    90    92    94    96]
%freq. 81Hz,83,85,87,89,91,93,95Hz
TXD = mean(mean(Tdr(binsM,:,:),3),1)';

binsR = binsM+1;
binsR = 1:100; 
binsR(binsM) = []; 
binsR(1:2) = []; 

FP = mean(mean(Tdr(binsR,:,:),3),1)';




%% mostrar resultados 
%clc
%1 - Taxa de detecção 

figure 
plot(TXD_semCorrecao(:,2),'.k','MarkerSize',6)
hold on 
plot([0 size(TXD_corrigido,1)],[TXD_corrigido(1,2) TXD_corrigido(1,2)], ':r','LineWidth',1.5)
%plot(TXD_corrigido(:,3),'*','Color',[0 0.4470 0.7410]','Markersize',6)
plot(TXD_corrigido(:,2),'.b','MarkerSize',6)
xlabel('Índice dos conjuntos de Parâmetros','fontsize',12)
ylabel('Taxa de Detecção(%)','fontsize',12)
box off

TXD_semCorrecao(:,3)=TXD;
TXD_corrigido(:,3)=TXD;


%2 - Falsos Positivo  
figure 
plot(FP_semCorrecao(:,2)*100,'.k','MarkerSize',8)
hold on 
plot(FP_corrigido(:,2)*100,'.b','Markersize',8)
%hold on 
plot([0 size(FP_corrigido,1)],[5 5], ':r','LineWidth',2)
xlabel('Índice dos conjuntos de Parâmetros','fontsize',12)
ylabel('Taxa de Falso Positivo(%)','fontsize',12)
box off
title('40 db dois resultados na mesma figura')

figure
plot(FP)
figure
boxplot(fpboxplot)
xlabel('NDC sem correção')

FP_semCorrecao(:,3)=FP;
FP_corrigido(:,3)=FP;

%começa aqui os subplots
figure
subplot(2,1,1) 
plot(FP_semCorrecao(:,1)*100,'.k','MarkerSize',5)
hold on 
plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
plot(FP_corrigido(:,1)*100,'.b','Markersize',5)
xlabel('Índice dos conjuntos de Parâmetros','fontsize',12)
ylabel('Taxa de Falso Positivo(%)','fontsize',10)
box off
title('Intensidade de 50 dB SPL','Fontsize',12)

subplot(2,1,1)
plot(FP_semCorrecao(:,3)*100,'.k','MarkerSize',5)
hold on 
plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
plot(FP_corrigido(:,3)*100,'.b','Markersize',5)
xlabel('Índice dos conjuntos de Parâmetros','fontsize',12)
ylabel('Taxa de Falso Positivo(%)','fontsize',10)
box off
title('Intensidade de 30 dB SPL','Fontsize',12)

subplot(3,1,3)
plot(FP_semCorrecao(:,3)*100,'.k','MarkerSize',7)
hold on 
plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
plot(FP_corrigido(:,3)*100,'.b','Markersize',7)
xlabel('Índice dos conjuntos de Parâmetros','fontsize',12)
ylabel('Taxa de Falso Positivo(%)','fontsize',10)
box off
title('30 db')

%termina aqui os subplots
figure
subplot(3,1,1) 
plot(FP_corrigido(:,1)*100,'.k','MarkerSize',9)
hold on 
plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
xlabel('Índice dos conjuntos de Parâmetros','fontsize',12)
ylabel('Taxa de Falso Positivo(%)','fontsize',10)
box off

subplot(3,1,2)
plot(FP_corrigido(:,2)*100,'.k','MarkerSize',9)
hold on 
plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
xlabel('Índice dos conjuntos de Parâmetros','fontsize',12)
ylabel('Taxa de Falso Positivo(%)','fontsize',10)
box off

subplot(3,1,3)
plot(FP_corrigido(:,3)*100,'.k','MarkerSize',9)
hold on 
plot([0 size(FP_semCorrecao,1)],[5 5], ':r','LineWidth',2)
xlabel('Índice dos conjuntos de Parâmetros','fontsize',12)
ylabel('Taxa de Falso Positivo(%)','fontsize',10)
box off


%%
% taxa de detecção x tempo 
timeM = time(binsM,:); 
timeM(timeM==-1) = Mmax;
timeM = mean(timeM,1)'*1; %1segundo por janela
TXD = mean(mean(Tdr(binsM,:,:),3),1)' *100;

%plot da fronteira parero
figure1 = figure; axes1 = axes('Parent',figure1); hold(axes1,'on');  %estudar
plot([0 1]*100, [Mmax Mmax],'-.k','linewidth',1) %linha do mmax 240
plot([TXD(1) TXD(1)], [min(timeM) max(timeM)],'-.b','linewidth',1) %linha do teste unico

figure
for ii = 1:size(parametros,1)
    plot(TXD(ii),timeM(ii),'.k','Markersize',6,'DisplayName',[num2str(parametros(ii,1)) '-' num2str(parametros(ii,2))])
end

[ p, idxs] = paretoFront([TXD,(-timeM)] ); 
auxL = p(:,1)<0.5; 
p(auxL,:) = [];
idxs(auxL,:) = [];

[~,ind] = sort(p(:,1));
p = p(ind,:);
idxs = idxs(ind,:);
% p([3,4,7,9],:) = [];
% idxs([3,4,7,9],:) = [];

plot(TXD(idxs),timeM(idxs),'-or','Markersize',8,'linewidth',1.2) %plot da fronteira

% ylim([min(-p(:,2))*.9 max(-p(:,2))*1.1])
% xlim([min(p(:,1))*80 max(p(:,1))*104])

set(axes1,'XMinorTick','on');
set(axes1,'YMinorTick','on');
box(axes1,'off');
hold off
xlabel('Detection Rate (%)','fontsize',12); 
ylabel('Mean Exam Time (min)','fontsize',12);
fprintf('\n'); 

for ii = 1:size(idxs,1)

    [I] = find((TXD == TXD(idxs(ii))) & (timeM==timeM(idxs(ii))));
    
    fprintf('\nPD = %f Tempo = %f ',TXD(idxs(ii)),timeM(idxs(ii))); 
    fprintf(' NI = %d ',length(I)); 
     I = I(1); 
    
    for jj = I  
        
        fprintf(' - Buffer:%d, M_step:%d', parametros(jj,1),parametros(jj,2)); 
    %    text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ','  num2str(parametros(jj,3)) '\}' ]);
%            text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ',100,'  num2str(parametros(jj,3)) '\}' ]);
%        text(PD_T(jj)*99,T_exame_G(idxs(ii))*.985,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) ',100\}' ]);
       text(TXD(jj),timeM(idxs(ii))*.975,['\{' num2str(parametros(jj,1)) ',' num2str(parametros(jj,2)) '\}' ]);
    end    
end
fprintf('\n'); 
xlim([min(TXD(idxs))*.95,max(TXD(idxs))*1.05])
ylim([min(timeM(idxs))*.95 Mmax*1.05])


