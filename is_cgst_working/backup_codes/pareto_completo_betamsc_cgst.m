% OBJ: gerar curva eficiencia de pareto para CGST MSC
% https://pt.wikipedia.org/wiki/Efici%C3%AAncia_de_Pareto

%% Setup
clearvars; close all; clc

%% Parametros de deteccao
FPd = 0.05;     % taxa de falso positivo desejado para o exame

% 1. Felix, Leonardo Bonato, et al.  (2006)
% "Statistical aspects concerning signal coherence applied to randomly 
% modulated periodic signals." IEEE Signal Processing Letters.
% 2. colatina
vec_volunt = 1:11; %4:7; %[1 2 3]; % [1 5 11]; % 1:11; % 1:5; %
% vec_intens = [1 3];  % [1 3 6]; %[1 3:6]; %1:6; % 0:5;  % 3:4; %
%                        % pulando intensidade 2 = 60 dB pois Mlimite = 50

% {'70dB';'60dB';'50dB';'40dB';'30dB';'ESP'}
vec_intens = [1 3 5]; % selecionando mesmos que colatina faster 2020


% Add citacao colatina
signal_freq_bins = [82    84    86    88    90    92    94    96];
noise_freq_bins  = 351:12:439; %351:1:451; 
freq_bins = [signal_freq_bins, noise_freq_bins];

%% Parametros de duracao:
% Em duvida, ver:
% "estudo_combinacoes_razoaveis_dados_colatina.m", e % "how_much_data.m"
Mlimite = [115,   58,  177,  298,  475,  477 ];
Mmin = 10; % 12 
Mmax = 58;
Mexpl= 10;
vec_M = [10 16 20 24 30 40]; %10:4:48; % Mmin:Mexpl:Mmax;
% vec_M = [10, 20, 40];

% limita para Mmax onde Mlimite>Mmax selecionado
Mlimite(Mlimite >Mmax) = Mmax;

%
Mlimite(Mlimite <min([vec_M,Mmin])) = Mmin;

Klimite_desejado = 50;
Kexpl = 5;
Kmin = [3,     3    ,   3,    3,    3,    3]; 
Kmax =  min(Mlimite,Klimite_desejado);
Mstep_min = floor(Mlimite./Kmax); % 1 sempre
Mstep_max = floor(Mlimite./Kmin); % depende to N max de teste e t max exame 

% Rever como padronizar o tamanho do vec_K, caso contrario as tabelas 
% terao tamanhos variaveis (quebra o codigo)
% entao, atualmente usando apenas intensidades que respeitam Kmax = 100
vec_K = [5,15,25,35,40];% Kmin(1):Kexpl:Kmax(1);
% Kmax = Klimite_desejado;
Kmax = max(vec_K);

disjoint_percentage = [0,10,20,30,40,50,80,100]./100;

Mstep_values = zeros(length(Mstep_max), length(disjoint_percentage));
for idx = 1:numel(disjoint_percentage)
    current_disjoint_percentage = disjoint_percentage(idx);
    switch current_disjoint_percentage
        case 1
            Mstep_values(:,idx) = Mstep_max;
        case 0
            Mstep_values(:,idx) = ones(size(Mstep_values(1)));
        otherwise
            Mstep_values(:,idx) = round(Mstep_max.*current_disjoint_percentage);
    end
end


%% Executar testes

t_decisao = zeros(...
    numel(vec_intens),...
    numel(vec_volunt),...
    numel(vec_M),...
    numel(Mstep_values),...
    Kmax,...
    numel(freq_bins)...
    );

TP = zeros(...
    numel(vec_intens),...
    numel(vec_volunt),...
    numel(vec_M),...
    numel(Mstep_values),...
    Kmax,...
    numel(freq_bins) ...
    );
FP = TP;
TN = TP;
FN = TP;

lambda = 0.3372*1.2;

disp(['Duracao minima estimada em minutos = '...
    num2str(numel(vec_intens)*numel(vec_volunt)*numel(vec_M)*...
    numel(Mstep_values)*numel(vec_K)*lambda/(60))])

disp(['Duracao minima estimada em horas = '...
    num2str(numel(vec_intens)*numel(vec_volunt)*numel(vec_M)*...
    numel(Mstep_values)*numel(vec_K)*lambda/(60*60))])

disp('_____________________________________')

disp(['Duracao pessimista estimada em minutos = '...
    num2str(numel(vec_intens)*numel(vec_volunt)*numel(vec_M)*...
    numel(Mstep_values)*numel(vec_K)*(lambda+lambda^2)/(60))])

disp(['Duracao pessimista estimada em horas = '...
    num2str(numel(vec_intens)*numel(vec_volunt)*numel(vec_M)*...
    numel(Mstep_values)*numel(vec_K)*(lambda+lambda^2)/(60*60))])

disp('_____________________________________')

disp(['Duracao P>95% estimada em minutos = '...
    num2str(numel(vec_intens)*numel(vec_volunt)*numel(vec_M)*...
    numel(Mstep_values)*numel(vec_K)*1.456/(60))])

disp(['Duracao P>95% estimada em horas = '...
    num2str(numel(vec_intens)*numel(vec_volunt)*numel(vec_M)*...
    numel(Mstep_values)*numel(vec_K)*1.456/(60*60))])

% PESSIMO!
% [t_decisao,TP,FP,TN,FN] = f_pareto(intens, volunt, M, %disj, K,freq_bins)
% Nao estou convencido que vale a pena 

tempos = nan(numel(vec_intens),numel(vec_volunt),numel(vec_M),...
    numel(Mstep_values),numel(vec_K));
t1 = tic;
for idx_intens = 1:numel(vec_intens)
    % approx 20seg max por loop
    for idx_volunt = 1:numel(vec_volunt)

        % Carrega dados (deveria ser programacao dinamica, custo atoa aqui)
        SIGNALS = terrible_load(vec_intens(idx_intens), ...
                        vec_volunt(idx_volunt)); %(intens, volunt);
        
        for idx_M = 1:numel(vec_M)
            for idx_Mstep = 1:numel(Mstep_values)
                for idx_K = 1:numel(vec_K) 
                    t2 = tic;

                    % Run CGST test
                    [t_decisaoi,TPi,FPi,TNi,FNi] = run_betamsc_cgst_RATE(...
                        vec_intens(idx_intens), ...
                        vec_volunt(idx_volunt), ...
                        vec_M(idx_M), ...
                        Mstep_values(idx_Mstep), ...
                        vec_K(idx_K), ...
                        Kmax,...
                        signal_freq_bins, ...
                        noise_freq_bins, ...
                        SIGNALS);
    
                    % Store results 
                     t_decisao(idx_intens, idx_volunt, idx_M, idx_Mstep,:,:) = t_decisaoi;
                     TP(idx_intens, idx_volunt, idx_M, idx_Mstep,:,:) = TPi;
                     FP(idx_intens, idx_volunt, idx_M, idx_Mstep,:,:) = FPi;
                     TN(idx_intens, idx_volunt, idx_M, idx_Mstep,:,:) = TNi;
                     FN(idx_intens, idx_volunt, idx_M, idx_Mstep,:,:) = FNi;

                tempos(idx_intens, idx_volunt, idx_M, idx_Mstep,idx_K) = toc(t2);
                end
            end
        end

        figure(1)
        histogram(tempos(:),'normalization','probability')
        disp(quantile(tempos(:),0.95))
        drawnow
        
    end
end
disp('Tempo total')
toc(t1) 

disp('Tempo por iteracao (exame)')
mean(tempos, 'all', 'omitnan')
%% Sanity check


% TP = zeros(...
%     numel(vec_intens),...
%     numel(vec_volunt),...
%     numel(vec_M),...
%     numel(Mstep_values),...
%     Kmax,...
%     numel(freq_bins) ...
%     );

TXD = sum_reduce(sum_reduce(sum_reduce(TP,6),1),1);
TXD = 100*TXD./(numel(vec_volunt)*numel(vec_intens)*numel(signal_freq_bins));


TXD_on_exam = mean_reduce(TXD, 3); % sum mean detection over diferent stages

TXD_mean_over_M = mean_reduce(TXD, 1); % sum mean detection over diferent stages

TXD_mean_over_all = mean_reduce(mean_reduce(mean_reduce(TXD, 1), 1)', 1);


%% Mesmo para fp

TXF = sum_reduce(sum_reduce(sum_reduce(FP,6),1),1);
TXF = 100*TXF./(numel(vec_volunt)*numel(vec_intens)*numel(noise_freq_bins   )) ;


TXF_on_exam = mean_reduce(TXF, 3); % sum mean detection over diferent stages

TXF_mean_over_M = mean_reduce(TXF, 1); % sum mean detection over diferent stages

TXF_mean_over_all = mean_reduce(mean_reduce(mean_reduce(TXF, 1), 1)', 1);


%% Para duracao

t = t_decisao;
t(t<0) = -1*t(t<0);
% t(t>0) = find(t>0);
% t(t<0) = find(t<0);
t(t==0) = NaN;

t = mean_reduce(mean_reduce(mean_reduce(t,6),1),1);
% t = 100*t./(numel(vec_volunt)*numel(vec_intens)*numel(signal_freq_bins)) 


t_on_exam = mean_reduce(t, 3); % sum mean detection over diferent stages

t_mean_over_M = mean_reduce(t, 2)'; % sum mean detection over diferent stages

t_mean_over_all = mean_reduce(mean_reduce(mean_reduce(t, 1), 1)', 1);
    

%%
figure
subplot(121)
plot(TXD(:), t(:),'b.', 'MarkerSize',10)
xlabel('Mean detection rate [%]', 'FontSize',16)
ylabel('Mean exam duration [s]', 'FontSize',16)
title('TP Pareto', 'FontSize',14)
grid on 

subplot(122)
plot(TXF(:), t(:),'rx', 'MarkerSize',10)
xlabel('Mean false positive rate [%]', 'FontSize',16)
ylabel('Mean exam duration [s]', 'FontSize',16)
title('FP Pareto', 'FontSize',14)
grid on
%%
figure
plot3(TXD(:),TXF(:),t(:),'.')
xlabel('Mean detection rate [%]', 'FontSize',16)
ylabel('Mean false positive rate [%]', 'FontSize',16)
zlabel('Mean exam duration [s]', 'FontSize',16)
grid on
%%
% TXF_mstep = mean_reduce(mean_reduce(TXF,3),1);
TXD_mstep = mean_reduce(mean_reduce(TXD,3),1);
t_mstep = mean_reduce(mean_reduce(t,3),1);
% plot3(TXF_mstep(:),Mstep_values(:),t_mstep(:),'.')
% [TXD_mstep(:),TXF_mstep(:),Mstep_values(:),t_mstep(:)]
% grid on

figure
subplot(121)
plot3(TXD_mstep(:),Mstep_values(:),t_mstep(:),'b.', 'MarkerSize',10)
xlabel('Mean detection rate [%]', 'FontSize',16)
ylabel('Number of samples per epoch [Mmin]', 'FontSize',16)
zlabel('Mean exam duration [mins]', 'FontSize',16)
title('TP Pareto', 'FontSize',14)
grid on

subplot(122)
plot3(TXF_mstep(:),Mstep_values(:),t_mstep(:),'rx', 'MarkerSize',10)
xlabel('Mean false positive rate [%]', 'FontSize',16)
ylabel('Sample/time increment  between tests [Mstep]', 'FontSize',16)
zlabel('Mean exam duration [mins]', 'FontSize',16)
title('FP Pareto', 'FontSize',14)
grid on

disp('TXDetec, TXFalsoP, Mstep, tempo_medio')
[TXD_mstep(:),TXF_mstep(:),Mstep_values(:),t_mstep(:)]

%%


TXF_m = mean_reduce(mean_reduce(TXF,2),2);
TXD_m = mean_reduce(mean_reduce(TXD,2),2);
t_m = mean_reduce(mean_reduce(t,2),2);

figure
subplot(121)
plot3(TXD_m(:),vec_M(:),t_m(:),'b.', 'MarkerSize',10)
xlabel('Mean detection rate [%]', 'FontSize',16)
ylabel('Number of samples per epoch [Mmin]', 'FontSize',16)
zlabel('Mean exam duration [mins]', 'FontSize',16)
title('TP Pareto', 'FontSize',14)
grid on

subplot(122)
plot3(TXF_m(:),vec_M(:),t_m(:),'rx', 'MarkerSize',10)
xlabel('Mean false positive rate [%]', 'FontSize',16)
ylabel('Number of samples per epoch [Mmin]', 'FontSize',16)
zlabel('Mean exam duration [mins]', 'FontSize',16)
title('FP Pareto', 'FontSize',14)
grid on

disp('TXDetec, TXFalsoP, Mmin, tempo_medio')
[TXD_m(:),TXF_m(:),vec_M(:),t_m(:)]

%%
TXF_k = mean_reduce(mean_reduce(TXF,1),1);
TXD_k = mean_reduce(mean_reduce(TXD,1),1);
t_k = mean_reduce(mean_reduce(t,1),1);

figure
subplot(121)
plot3(TXD_k(:),1:Kmax,t_k(:),'b.', 'MarkerSize',10)
xlabel('Mean detection rate [%]', 'FontSize',16)
ylabel('Number of tests [K]', 'FontSize',16)
zlabel('Mean exam duration [mins]', 'FontSize',16)
title('TP Pareto', 'FontSize',14)
grid on

subplot(122)
plot3(TXF_k(:),1:Kmax,t_k(:),'rx', 'MarkerSize',10)
xlabel('Mean false positive rate [%]', 'FontSize',16)
ylabel('Number of tests [K]', 'FontSize',16)
zlabel('Mean exam duration [mins]', 'FontSize',16)
title('FP Pareto', 'FontSize',14)
grid on

disp('TXDetec, TXFalsoP, K, tempo_medio')
[TXD_k(:),TXF_k(:),[1:Kmax]',t_k(:)]


%% Salvar resultados
save('backup5_pareto_completo_betamsc_cgst.mat')

%%
% Get the linear index of the 
linearIndex = find(TXF>5);
% Convert linear index to subscript
[b_M, b_Mstep, b_Kmax,vals] = ind2sub(size(TXF), linearIndex);
figure
subplot(131)
histogram(b_Mstep, 'normalization', 'probability')
subplot(132)
histogram(b_Kmax, 'normalization', 'probability')
subplot(133)
histogram(b_M, 'normalization', 'probability')

%%

t_pareto = t_decisao;
t_pareto(t_pareto<0) = -1*t_pareto(t_pareto<0);
t_pareto(t_pareto==0) = NaN;

t_pareto = mean_reduce(mean_reduce(t_pareto,6),2);
TXD_pareto = sum_reduce(sum_reduce(TP,6),2);
TXD_pareto = 100*TXD_pareto./(numel(vec_volunt)*numel(signal_freq_bins));

TXD_pareto = TXD_pareto(2,:,:,:);
t_pareto = t_pareto(2,:,:,:);

[ p, idxs] = paretoFront([TXD_pareto(:),(-t_pareto(:))] ); 
auxL = p(:,1)<0.5; 
p(auxL,:) = [];
idxs(auxL,:) = [];
[~,ind] = sort(p(:,1));
p = p(ind,:);
idxs = idxs(ind,:);
% p([3,4,7,9],:) = [];
% idxs([3,4,7,9],:) = [];
idxs_2006 = idxs;
x_plot = TXD_pareto(idxs_2006);
y_plot = t_pareto(idxs_2006);
plot(y_plot(:),x_plot(:),'c-o','Markersize',8,'linewidth',1.8) 

%%
load('C:\PPGEE\Assessing CGST on ASSR\new_code\backup3_pareto_completo_betamsc_cgst.mat')
%%


t_pareto = t_decisao;
t_pareto(t_pareto<0) = -1*t_pareto(t_pareto<0);
t_pareto(t_pareto==0) = NaN;

t_pareto = mean_reduce(mean_reduce(t_pareto,6),2);
TXD_pareto = sum_reduce(sum_reduce(TP,6),2);
TXD_pareto = 100*TXD_pareto./(numel(vec_volunt)*numel(signal_freq_bins));

TXD_pareto = TXD_pareto(2,:,:,:);
t_pareto = t_pareto(2,:,:,:);
% 
% TXD_pareto(isnan(TXD_pareto)) = [];
% t_pareto(isnan(t_pareto)) = [];

[ p, idxs] = paretoFront([TXD_pareto(:),(t_pareto(:))] ); 
auxL = p(:,1)<0.05; 
p(auxL,:) = [];
idxs(auxL,:) = [];

nan_rows = isnan(p(:,2));
row_indices = find(nan_rows);
p = p(~isnan(p(:,2)), :);

[~,ind] = sort(p(:,1));
p = p(ind,:);
idxs = idxs(ind,:);

x_plot = TXD_pareto(idxs);
y_plot = t_pareto(idxs);

% uiopen('C:\PPGEE\Assessing CGST on ASSR\Pareto Todos.fig',1)
%uiopen('C:\PPGEE\Assessing CGST on ASSR\Pareto Todos_30.fig',1)
uiopen('C:\PPGEE\Assessing CGST on ASSR\Pareto Todos_50.fig',1)
%uiopen('C:\PPGEE\Assessing CGST on ASSR\Pareto Todos_70.fig',1)

hold on
plot(y_plot(:),x_plot(:),'c-o','Markersize',8,'linewidth',1.8) 




a=gca;
a.Legend.String = {'(Zonateli, 2020)',  '(Bazoni,2021)', '(Antunes, 2019)',...
    '(Vaz, 2023)','Single Test','Beta-CGST'};
a.Legend.AutoUpdate = 'off';

scatter(t_pareto(:), TXD_pareto(:),8, 'filled','MarkerEdgeColor','cyan',...
              'MarkerFaceColor','cyan')

% ylim([0,480])
% xlim([0,40])
ylim([0,255])
xlim([0,75])

%%
sum(TXF<5)/numel(TXF)
sum(TXF<5, 'all')/numel(TXF)
whos TXF
TXF2 = sum_reduce(sum_reduce(FP,6),2);
whos TXF2
TXF2_50dB = TXF2(2,:,:,:);
sum(TXF2_50dB<5, 'all')/numel(TXF2_50dB)
TXF2_30dB = TXF2(3,:,:,:);
sum(TXF2_30dB<5, 'all')/numel(TXF2_30dB)
TXF2_70dB = TXF2(1,:,:,:);
sum(TXF2_70dB<5, 'all')/numel(TXF2_70dB)

%%
% intensIv = []; tamMv = []; tamMstepv = []; numKv = [];
[intensIv, tamMv, tamMstepv, numKv] = ind2sub(size(TXF2),find(TXF2>5));


figure(1)
subplot(2,2,1)
histogram(intensIv,size(TXF2,1),'normalization','cdf')
xlabel('30 50 70 dB')
title('intensIv')

subplot(2,2,2)
histogram(tamMv,size(TXF2,2)+1,'normalization','cdf', 'BinEdges',1:numel(vec_M))
xlabel('10    16    20    24    30')
title('tamMv')

subplot(2,2,3)
histogram(tamMstepv, size(TXF2,3),'normalization','cdf')
title('tamMstepv')

subplot(2,2,4)
histogram(numKv, size(TXF2,4),'normalization','cdf','BinEdges', vec_K)
title('numKv')


figure(2)
subplot(2,2,1)
histogram(intensIv,size(TXF2,1),'normalization','probability')
xlabel('30 50 70 dB')
title('intensIv')

subplot(2,2,2)
histogram(tamMv,size(TXF2,2),'normalization','probability', ...
    'BinEdges',1:numel(vec_M))
xlabel('10    16    20    24    30    40')
title('tamMv')

subplot(2,2,3)
histogram(tamMstepv, size(TXF2,3),'normalization','probability')
title('tamMstepv')

subplot(2,2,4)
histogram(numKv, size(TXF2,4),'normalization','probability', ...
    'BinEdges', vec_K)
title('numKv')
