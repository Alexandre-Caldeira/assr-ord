% OBJ: gerar curva eficiencia de pareto para CGST MSC
% https://pt.wikipedia.org/wiki/Efici%C3%AAncia_de_Pareto

%% Setup
clearvars; close all; clc

%% Parametros
Kmax = 15;      % numero maximo de testes aplicados em sequencia 
FPd = 0.05;     % taxa de falso positivo desejado para o exame

% Felix, Leonardo Bonato, et al.  (2006)
% "Statistical aspects concerning signal coherence applied to randomly 
% modulated periodic signals." IEEE Signal Processing Letters.

vec_volunt = 1:11; % 1:11; % 1:5; %
vec_intens = [0:5]; % 1:5;  % 3:4; %

%%%%%%%%%%%%
signal_freq_bins =  [82   90    84    86    88    90    92    94    96];
noise_freq_bins = 5; % round(signal_freq_bins.*exp(1)/2)+5
%%%%%%%%%%%%


% 
Mlimite = nan(length(vec_intens), length(vec_volunt), 2);
for intens = vec_intens
    for volunt = vec_volunt
        %% Carrega dados
        
        % SET THIS PATH:
        path = 'C:\PPGEE\SBEB_CBA_24\ASSR - Coleta OFFLINE';
        addpath(path)
        %vetor dos voluntários 
        Vvoluntario = {'Abdon';'Ana';'BBB';'Colatina';'Erick';'Luciana';...
            'Sombra';'Quenaz';'Vinicius';'Sacola';'Wreikson'}; 
        
        %vetor da intensidade 
        Vintensidade = {'70';'60';'50';'40';'30'}; 
        load('eletrodos.mat')
        ganho  = 200;
        remoc = 0.1/ganho; 
        
        voluntario = cell2mat(Vvoluntario(volunt,:));
        
        if intens==0
             load([voluntario 'ESP'], 'x','Fs','binsM','freqEstim') 
        else
            intensidade = cell2mat(Vintensidade(intens,:));
            load([voluntario '_'  intensidade 'dB'], 'x','Fs','binsM','freqEstim') 
        end
        
        nfft = Fs;
        x = x - repmat(mean(x),nfft,1);
        x(:,1:2,:) =[]; 
        
        % %encontrar o valor máximo por canal 
        Vmax = squeeze(max(max(abs(x)),[],3));
        ind = Vmax>remoc;
        
        pos_eletrodo= 1;
        xmedia = x(:,~ind,pos_eletrodo);
        
        SIGNALS = fft(xmedia,Fs);%*2/nfft*1e9;
        FS = Fs;
        NFFT = Fs;
        SIGNALS = SIGNALS(1:floor(end/2)+1,:); 
        Y = SIGNALS;
        Y(59:63,:) = 0;
        all_freqs = [signal_freq_bins noise_freq_bins];

        Mlimite(intens+1, volunt, (Fs==1000)+1) = size(SIGNALS,2);

    end
end

%%
disp(' Analisando percentis: ')
percentis = [.2,.25,.3,.4,.6,.8,.95];
fprintf('%d\t',round(100*percentis))
disp(' ')

disp('_ Para Fs = 1000 Hz _')
Mstat = nan(length(vec_intens), length(vec_volunt));
Mstat = Mlimite(:,:,1);

% duracao_minima_por_intensidade = min(Mstat,[],2,"omitnan")
% duracao_media_por_intensidade = mean(Mstat,2,"omitnan")
% duracao_maxima_por_intensidade = max(Mstat,[],2,"omitnan")

duracao_pct_por_intensidade_1000hz = quantile(Mstat,percentis,2)
desvio_padrao_1000Hz = std(duracao_pct_por_intensidade_1000hz,0,2)

disp('_ Para Fs = 1750 Hz _')
Mstat = nan(length(vec_intens), length(vec_volunt));
Mstat = Mlimite(:,:,2);

% duracao_minima_por_intensidade = min(Mstat,[],2,"omitnan")
% duracao_media_por_intensidade = mean(Mstat,2,"omitnan")
% duracao_maxima_por_intensidade = max(Mstat,[],2,"omitnan")

duracao_pct_por_intensidade_1750hz = quantile(Mstat,percentis,2)
desvio_padrao_1750Hz = std(duracao_pct_por_intensidade_1750hz,0,2)

disp('_ Para conjunto total de Fs _')
Mstat = nan(length(vec_intens), length(vec_volunt));
Mstat = Mlimite(:,:,1);
aux = Mlimite(:,:,2);
Mstat(isnan(Mstat)) = aux(isnan(Mstat));

% duracao_minima_por_intensidade = min(Mstat,[],2)
% duracao_media_por_intensidade = mean(Mstat,2)
% duracao_maxima_por_intensidade = max(Mstat,[],2)

duracao_pct_por_intensidade_todos = quantile(Mstat,percentis,2)
desvio_padrao_todos = std(duracao_pct_por_intensidade_todos,0,2)


% Escolha de Mlimite pode ser feita como menor valor que ultrapassa 1 sigma
num_maiores_que_desvio_padrao = sum(...
    duracao_pct_por_intensidade_todos-sqrt(var(duracao_pct_por_intensidade_todos,0,2))>0,1);

Mlimite_por_intensidade= floor(duracao_pct_por_intensidade_todos(:,...
    find(num_maiores_que_desvio_padrao==length(vec_intens),1,'first')))'
Percentil_contemplado = 100*(1-...
    percentis(find(num_maiores_que_desvio_padrao==length(vec_intens),1,'first')))

Mlimite_por_intensidade =   [115    58   177   298   475   477];


%%
% 
% Mmin = 12;       
% possible_Mstep = [1,2,4,6,8,10,12,16,24,28,30,32,36,40];
% Mmax = 64;
% 
% for current_Mstep = possible_Mstep
%     % disp(size(Mmin:current_Mstep:Mmax))
%     disp(floor((Mmax-Mmin)/current_Mstep))
% end
