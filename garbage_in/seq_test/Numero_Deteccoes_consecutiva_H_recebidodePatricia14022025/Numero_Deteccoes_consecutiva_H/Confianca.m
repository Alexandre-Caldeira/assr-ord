% Intervalo de Confiança
clear all

%NN = 1000; % nº amostras
vol = 11;
n1 = 8; % qtidade de frequencias
n2 = (104 -70-8+1);
N1 = vol * n1; %Espontâneo
N2 = vol * n2; %Bandas Laterais
%% 
p = 0.05; %nível de significãncias 
%limite_conf_sup = 1- p;
limite_conf_sup = .90; %intervalo de confiança
%se aumentar para 0.95 -> para de 90% limite de 92%
%se aumentar para 0.90 -> limite em simulação de 85%.
pp = (1-limite_conf_sup)/2; 

%NN = N1;
NN = N2;
Nrun = 100000;
for ii = 1:Nrun
    r(ii) = binornd(NN,p); %obrigatorio nivel de signific}anica
end

ni = quantile(r,1-limite_conf_sup);
ns = quantile(r,limite_conf_sup);
mean(r>=ni);
mean(r<=ns);
mean((r<=ns).*(r>=ni));


% linf = (quantile(r,1-limite_conf_sup))/N1;
% lsup = (quantile(r,limite_conf_sup))/NN;
% lsup + linf;

linf = (quantile(r,pp)+1)/NN;
lsup = (quantile(r,1-pp)-1)/NN;

n_sinais = [1 2 4 6 8];

margem_erro = 100*lsup/NN;

% ni = quantile(r,1-limite_conf_sup) %nº de amostras inferior [discreto]
% ns = quantile(r,limite_conf_sup) %nº de amostras superior [discreto]
% mean(r>=ni)
% mean(r<=ns)
% mean((r<=ns).*(r>=ni)) % Intervao de confiança 
% linf = (quantile(r,1-limite_conf_sup))/N1 %Limite inferior
% lsup = (quantile(r,limite_conf_sup))/N1 %Limite superior
% % 
% lsup + linf;
%% 
% %limite_conf_sup = 0.6; %0,6
%     pp = (1-limite_conf_sup)/2; 
%     
%     ni = (quantile(r,pp));
%     ns = (quantile(r,1-pp));
%     mean(r>=ni);
%     mean(r<=ns);
%     mean((r<=ns).*(r>=ni))
% 
% %    linf = (quantile(r,pp)+1)/NN;
%      linf = (quantile(r,pp))/N1;
%  %   lsup = (quantile(r,1-pp)-1)/NN;  
%     lsup = (quantile(r,1-pp))/N1; 
