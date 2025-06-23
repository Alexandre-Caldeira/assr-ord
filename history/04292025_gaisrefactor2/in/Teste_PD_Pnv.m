%% 1.Gerar Curva PD
%************* Aqui TEM que alterar Mmax é 02 locais**********
%% 2. Determinar SNR para PD = 50%
clear all, clc
tic
Mmax = 20; % quantidades de janelas
L = 1000; % tamanho da janela
bin = 65; % frequência para evitar espalhamento
nRuns = 10000; %quantidade de testes
alfa = 0.05;
k = 1:L*Mmax; 

%valor teorico 
VC_teorico= 1 - alfa^(1/(Mmax-1));% tem ponto depois de alfa?

%Gerar Sinal (senoide+ruido)
MSC = zeros(nRuns,1); % Inicializar
TD = zeros(length(-58:-34),1);
aux = 1; % incrementar em TD
for SNR_db = -58:-34;
    disp (SNR_db);
    SNR = 10^(SNR_db/10);
    A = sqrt(L*SNR);
    %s = A * sin(2*pi*(bin - 1)/L*k);
    for i = 1:nRuns
        % Valores da DFT
        Y = A +(randn(Mmax,1) + 1j*randn(Mmax,1));
        % Cálculo da MSC
        msc =  abs(sum(Y)).^2./(Mmax*sum(abs(Y).^2));
        %Comparar com VC
        MSC(i) = msc > VC_teorico; % vetor de 0 e 1 da resposta da msc
    end
    %Taxa de detecção
    TD(aux) = mean (MSC);
    aux = aux + 1;
end

% Plotar PD
plot (-58:-34, TD);
xlabel('SNR','fontsize',12)
ylabel('Probabilidade de Detecção','fontsize',12)
grid on

% Determinar SNR para PD = 50%
Nii = 100;
msc = ones(Mmax,nRuns);
delta_SNR = 1;
passo = 0.75;
TD_ant = 0;
SNR_db = -43; 
for ii = 1:Nii
        
    SNR = 10^(SNR_db/10);
    A = sqrt(L*SNR);
    %s = A * sin(2*pi*(bin - 1)/L*k);
    
    for i = 1:nRuns
        % Valores da DFT
        Y = A +(randn(Mmax,1) + 1j*randn(Mmax,1));
        % Cálculo da MSC
        msc =  abs(sum(Y)).^2./(Mmax*sum(abs(Y).^2));        
        %Comparar com VC 
        MSC(i) = msc > VC_teorico; % vetor de 0 e 1 da resposta da msc
    end
    %Taxa de detecção
    TD = mean (MSC);
   
    if TD == 0.5
        break
    elseif TD > 0.5 && TD_ant > 0.5
        SNR_db = SNR_db - delta_SNR;
    elseif TD > 0.5 && TD_ant < 0.5
        delta_SNR = passo * delta_SNR;
        SNR_db = SNR_db - delta_SNR;
    elseif TD < 0.5 && TD_ant < 0.5
        SNR_db = SNR_db + delta_SNR;
    else TD < 0.5 && TD_ant > 0.5;
        delta_SNR = passo * delta_SNR;
        SNR_db = SNR_db + delta_SNR;
    end
    TD_ant = TD;
    disp([num2str(ii*100/Nii),'%'])
end
save(['SNR_Mmax_' num2str(Mmax) '_alfa_'  num2str(alfa) '_Tjanela_'  num2str(L) '.mat'],'SNR_db') 

%% MSC = SNR -42.9426 dB
clear all
clc
tic

L = 1000; % tamanho da janela
Mmax = 20; % quantidades de janelas
bin = 65; % frequência para evitar espalhamento
nRuns = 10000; %quantidade de testes 
k = 1:L*Mmax;
alfa = 0.05;
alfa_desejado = 0.05;
% VC_teorico= 1 - alfa^(1/(Mmax-1)); %valor teorico
load(['SNR_Mmax_' num2str(Mmax) '_alfa_'  num2str(alfa) '_Tjanela_'  num2str(L) '.mat'],'SNR_db')

% SNR_db = -42.9426; % Depois salvar e fazer load
disp (SNR_db);
SNR = 10^(SNR_db/10);
A = sqrt(L*SNR);
%%
%Teste sequencial  sinal+ruido
for i = 1:nRuns    
    Y = A +(randn(Mmax,1) + 1j*randn(Mmax,1));
    % Cálculo da MSC de forma sequencial
    NUM = cumsum(Y);
    DEN = (1:Mmax)'.*cumsum(abs(Y).^2);
    msc(:,i) = (abs(NUM).^2)./DEN;
end
% Verificar para CADA janela quais msc > VC_detecção
% alfa_teste = 0.001

M_2 = (2:Mmax)';
% buffer = 30;
%%
% Carregando o alfa inicial/teste encontrado em Teste_Seq_Pnv
% Esse alfa de detection para encontrar o VC_not
load(['Alfa_Corrigido_Detection_Mmax_' num2str(Mmax) '_alfa_'  num2str(alfa_desejado)  '_Tjanela_' num2str(L) '.mat'],'alfa')
%alfa = 0.001;  %alfa de cada teste 

VC_teste = 1-alfa.^(1./(M_2-1));
VC_teste = [1;VC_teste]; %acrescetou Um valor em VC_teste
%%

VC_not = ones (Mmax,1);
for i = 2:Mmax
    %Comparar com VC
    msc_gt = any(bsxfun (@gt,msc (i:Mmax,:),VC_teste(i:Mmax)));
    msc_lt = msc (i,:) < VC_teste (i);
    VC_not(i) = quantile(msc(i,msc_gt & msc_lt),alfa); %Encontrar percentil 1% e 99%
    
end
alfa = 0.05;          
save(['VC_not_Unitario_Mmax_' num2str(Mmax) '_alfa_'  num2str(alfa) '.mat'],'VC_not') 
%%
plot(VC_teste);
xlabel('Valores Críticos','fontsize',12)
ylabel('Janelas','fontsize',12)
grid on
hold on
plot(VC_not);
toc