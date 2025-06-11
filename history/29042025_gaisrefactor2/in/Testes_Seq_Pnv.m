clear all, clc

% 1.Curva de Valores críticos
M_2 = 2:160;
alfa = 0.01;
VC = 1-alfa.^(1./(M_2-1));

lw  =2; 
plot(M_2,VC,'r','LineWidth',lw) 
xlabel('Janelas','fontsize',12)
ylabel('Valores Críticos','fontsize',12)
grid on 


%% 2.Estimar valor crítico de MC
clear all
clc

L = 1024; % tamanho da janela
M_2 = 160; % quantidades de janelas
bin = 65; % frequência para evitar espalhamento
nRuns = 100; %quantidade de testes
alfa = 0.05;

%valor teorico 
VC_teorico= 1 - alfa^(1/(M_2-1));
msc = ones(M_2,nRuns);
for i = 1:nRuns
    y = randn(L*M_2,1); % sinal/ruido
    %Dividir em janelas (sinal,linhas,colunas) Sinal janelado
    yj = reshape(y,L,M_2);
    %Aplicar fft
    Y =fft(yj);
    
    for Mj = 2:1:M_2
        
        %Aplicar MSC, cada coluna vai ser um ruido
        msc (Mj,i) =  abs(sum(Y(bin,1:Mj),2)).^2./(Mj*sum(abs(Y(bin,1:Mj)).^2,2));
        
    end
    
end
plot (msc);

%% 3.Curva de Valores críticos
clear all
clc

L = 1024; % tamanho da janela
Mmax = 160; % quantidades de janelas
bin = 65; % frequência para evitar espalhamento
nRuns = 1000; %quantidade de testes
M_2 = (2:Mmax)';

alfa = 0.005;  %alfa teste 
alfa_desejado = 0.01;
VC_teste = 1-alfa.^(1./(M_2-1));
VC_teste = [1;VC_teste]; %acrescetou Um valor em VC_teste
buffer = 30;

% %valor teorico 
% VC_teorico= 1 - alfa^(1/(M-1));
msc = ones(Mmax,nRuns);
%MSC = ones(Mmax,nRuns);
for i = 1:nRuns
    y = randn(L*Mmax,1); % sinal/ruido
    %Dividir em janelas (sinal,linhas,colunas) Sinal janelado
    yj = reshape(y,L,Mmax);
    %Aplicar fft
    Y =fft(yj);
    
    for Mj = 2:1:Mmax
        
        %Aplicar MSC, cada coluna vai ser um ruido
        msc (Mj,i) =  abs(sum(Y(bin,2:Mj),2)).^2./(Mj*sum(abs(Y(bin,2:Mj)).^2,2));
        
    end
 
end
%Comparar com VC / any:verifica se tem UM
fp = any(bsxfun (@gt,msc (buffer:Mmax,:),VC_teste(buffer:Mmax)));
% for i = 1: nRuns
%    tp(i) = any (msc (buffer:Mmax,i) > VC_teste(buffer:Mmax));
% end
FP = mean(fp);
delta_alfa = alfa/4;
delta_crit = 0;
passo = 0.5; 

%% Ajuste de alfa_teste
clear all, clc
tic
L = 1000; % tamanho da janela
Mmax = 20; % quantidades de janelas
alfa_desejado = 0.05; 
bin = 65; % frequência para evitar espalhamento
nRuns = 10000; %quantidade de testes
M_2 = (2:Mmax)';
buffer =2; %Mmin

msc = ones(Mmax,nRuns);
for i = 1:nRuns
    Y = randn(Mmax,1) + 1j*randn(Mmax,1);
    % Cálculo da MSC
    NUM = cumsum(Y);
    DEN = (1:Mmax)'.*cumsum(abs(Y).^2);
    msc(:,i) = (abs(NUM).^2)./DEN;
end
alfa = 0.005;  %alfa teste 

Nii = 100;  % nº iterações para achar FP
delta_alfa = alfa/4;
passo = 0.75;
FP_ant = 0.05;
for ii = 1:Nii
        
    VC_teste = 1-alfa.^(1./(M_2-1));
    VC_teste = [1;VC_teste]; %acrescetou Um valor em VC_teste

    %Comparar com VC; any:verifica se tem UM
    fp = any(bsxfun (@gt,msc (buffer:Mmax,:),VC_teste(buffer:Mmax)));
    FP = mean(fp);
    
    if FP == alfa_desejado
        break
    elseif FP > alfa_desejado && FP_ant > alfa_desejado
        alfa = alfa - delta_alfa;
        if alfa < 0 
            alfa = 0;
        end
    elseif FP > alfa_desejado && FP_ant < alfa_desejado
        delta_alfa = passo * delta_alfa;
        alfa = alfa - delta_alfa;
        if alfa < 0
            alfa = 0;
        end
    elseif FP < alfa_desejado && FP_ant < alfa_desejado
        alfa = alfa + delta_alfa;
    else FP < alfa_desejado && FP_ant > alfa_desejado;
        delta_alfa = passo * delta_alfa;
        alfa = alfa + delta_alfa;
    end
    FP_ant = FP;
end
toc
% lw  = 2; %
%salvar as variáveis 
save(['Alfa_Corrigido_Detection_Mmax_' num2str(Mmax) '_alfa_'  num2str(alfa_desejado) '_Tjanela_' num2str(L) '.mat'],'alfa')

% plot(VC);
% xlabel('Valores Críticos','fontsize',12)
% ylabel('Janelas','fontsize',12)
% grid on
% hold on