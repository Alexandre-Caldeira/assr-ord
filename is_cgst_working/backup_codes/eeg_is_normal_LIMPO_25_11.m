%% Setup
clearvars; 
close all
clc

should_plot = 1;

%% Parametros
Kd = 999;        % numero total de testes a aplicar sequencialmente 
                % (dos ultimos Mmin segundos apos passar Msample segundos)
FPd = 0.05;     % taxa de falso positivo desejado para o exame

% Msample = 32; % best = 32 / 20 a 24
Mmin = 2; %floor(Msample*0.6);
sum_step = 1;

NumT = 5e3;
Nsinais = NumT;

int_inic = 3;
max_int = 3;

FNRt = nan(11,max_int-int_inic+1, Kd);
TPRt = nan(11,max_int-int_inic+1, Kd);
TNRt = nan(11,max_int-int_inic+1, Kd);
FPRt = nan(11,max_int-int_inic+1, Kd);

for cont_int = int_inic:max_int
for cont_vol=3:3

% always plot last run
if cont_int ==max_int
    if cont_vol==11
        should_plot = 1;
    end
end



    cont_vol

% REAL: s frequências das portadoras para ambas os ouvidos foram as mesmas: 
% 500, 1000, 2000,4000 Hz, modulados, respectivamente, nas frequências a
% 81, 85, 89 e 93 Hz, para o ouvido direito, e 83, 87, 91 e 95 Hz, 
% para o ouvido esquerdo . - pg 52 tese colatina

signal_freq_bins =  [82  84  86  88    90    92    94    96];
% noise_freq_bins = 300:1:300+Msample;%440:1:451;

% noise_freq_bins  =  [82  84  86  88    90    92    94    96];
% signal_freq_bins = 351:1:451;

% 443:1:451;
noise_freq_bins = 300+round(signal_freq_bins.*exp(1)/2)+5; %[300 400 500]; %
%%%%%%%%%%%%
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

voluntario = cell2mat(Vvoluntario(cont_vol,:));

if cont_int<0
     load([voluntario 'ESP'], 'x','Fs','binsM','freqEstim') 
else
    intensidade = cell2mat(Vintensidade(cont_int,:));
    load([voluntario '_'  intensidade 'dB'], 'x','Fs','binsM','freqEstim') 
end

nfft = Fs;%1segundo de sinal 
 
%retirar componente DC por janela (fiz isso pq no processamento em
%tempo real é por janela)
x = x - repmat(mean(x),nfft,1); %tirar a média de cada de cada trecho - devido a remoção

% meu
% fcInferior = 70;
% fcSuperior = 100;


% d1 = designfilt('bandstopiir','FilterOrder',4, ...
%                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
%                'DesignMethod','butter','SampleRate',Fs);


% Em seguida, os sinais foram filtrados digitalmente por um filtro
% Butterworth passa-banda de oitava ordem com frequência de corte de 1 Hz 
% acima e abaixo da frequência de modulação (Antunese Felix, 2019).
% - pg 53, tese colatina
% fcInferior = min(signal_freq_bins)-1;
% fcSuperior = max(signal_freq_bins)+1;
% d2 = designfilt('bandpassiir','FilterOrder',8, ...
%     'HalfPowerFrequency1',fcInferior,'HalfPowerFrequency2',fcSuperior, ...
%     'DesignMethod','butter','SampleRate',Fs);

% x = filtfilt(d1,x);
% x = filtfilt(d2,x);

%colatina
% fcInferior = 70; % 70
% fcSuperior = Fs/2 -1; % 100
% [b,a] = butter(8,[fcInferior/(Fs/2), fcSuperior/(Fs/2)]);
% x = filter(b,a,x); 
% excluir os dois primeiros segundos do inicio da coleta 
x(:,1:2,:) =[]; 

% %encontrar o valor máximo por canal 
Vmax = squeeze(max(max(abs(x)),[],3));
ind = Vmax>remoc;
% 
% xmedia = squeeze(mean(x(:,~ind,:),2));

  %   eletrodos =
  % 
  % 16×2 char array
  % 
  %   'FC'
  %   'F4'
  %   'T6'
  %   'P4'
  %   'T4'
  %   'Oz'
  %   'C4'
  %   'T5'
  %   'P3'
  %   'F7'
  %   'F3'
  %   'T3'
  %   'C3'
  %   'Fz'
  %   'Pz'
  %   'Cz'
% pos_eletrodo= 5;
pos_eletrodo= 1;
xmedia = x(:,~ind,pos_eletrodo);

SIGNALS = fft(xmedia,Fs);%*2/nfft*1e9;
FS = Fs;
NFFT = Fs;
SIGNALS = SIGNALS(1:floor(end/2)+1,:); % only half the FFT spectrum is valid
% f = FS/2*linspace(0,1,NFFT/2+1)'; % only half the FFT spectrum is valid
% max_length = size(SIGNALS,2);
Y = SIGNALS;
Y(59:63,:) = 0;
all_freqs = [signal_freq_bins noise_freq_bins];
% SIGNALS = normalize(SIGNALS(:,3:end));
% SIGNALS = normalize(SIGNALS);
% SIGNALS = SIGNALS./std(SIGNALS,0,2);

%% Calc thresholds
% M = floor(size(SIGNALS,2)/Kd);
if Kd>size(SIGNALS,2)
    if floor(size(SIGNALS,2)/Mmin)<=1
        % exame_tem_apenasfloor(size(SIGNALS,2)/Kd))

        % error
        disp(['Exame curto demais para teste sequencial com MSC! ', ...
           num2str(size(SIGNALS,2)),' janelas totais.' ])
        continue
    end  

    % K = floor(size(SIGNALS,2)/Mmin);
    K = size(SIGNALS,2);
    M=Mmin;
else
    K = Kd;
    M = Mmin;
    SIGNALS = SIGNALS(:,1:K);
end

Y = SIGNALS;
Mmin = 32;% 32; usually 32, n minimo de amostras para calc MSC. coincide com tamanho da janela de teste
NFFT = 1000;
Njanelas = K+1;
NpontosTotal = NFFT*Njanelas;

% disjoint:
% type = 'disjoint';
% Mstep = Mmin;
% Kmax_disj = floor(Njanelas/Mmin);
% Kmax = Kmax_disj;
% upperlim = Kmax -1;
% step = 1;

% half-repeating
type = '50% resample';
Mstep = round(Mmin/2);
Kmax_halfrep = round(Njanelas-Mmin/2);
Kmax = Kmax_halfrep;
upperlim =  numel(1:Mstep:Njanelas-Mmin);
step = 1;

% 1-out
% type = '1-out after min';
% Mstep = 1;
% Kmax_1step = Njanelas-Mmin; % = numel(tkf)-1 = numel(tki)
% Kmax = Kmax_1step;
% upperlim = Njanelas-Mmin;
% step = round((Njanelas-Mmin)/Kmax);


% Maxtests
% type = 'every window';
% Mstep = 1;
% Kmax = Njanelas;
% upperlim = Njanelas;
% step = round((Njanelas)/K);


tki = 1:Mstep:Njanelas-Mmin;
tkf = Mmin:Mstep:Njanelas;

%% resulting dist

for idx = 1:numel(tki)
    % figure(idx)
    % stem(abs(Y(:,tki(idx):tkf(idx))))
    
    rolling_msc(:,idx) = msc_fft(Y(:,tki(idx):tkf(idx)),Mmin);
    rolling_beta(:,idx) = betarnd(1,Mmin-1,1,NFFT);
end

CRMSC = cumsum(rolling_msc,2);
CRBETA = cumsum(rolling_beta,2);

mu_crmsc = zeros(1,numel(tki));
mu_crbeta = zeros(1,numel(tki));

normality_crmsc = zeros(1,numel(tki));
normality_crbeta = zeros(1,numel(tki));
    
for idx = 1:numel(tki)
    % a = normalize(CRMSC(:,idx));
    a = normalize(CRMSC(1:100,idx));
    % a = normalize(CRMSC(1:100,idx));
    % a(a>5*std(a)) = [];
    b = normalize(CRBETA(:,idx));

    mu_crmsc(idx) = mean(a);
    mu_crbeta(idx) = mean(b);

    figure(1)
    subplot(121)
    h1 = histogram( a , 100);
    normality_crmsc(idx) = 1 - kstest(a,'alpha', (5/100)/Kmax);

    title(['Test k = ', num2str(idx),' of ', num2str(numel(tki))])
    xlabel('CRMSC')
    ylabel('count')

    subplot(122)
    h2 = histogram( b , 100);
    normality_crbeta(idx) = 1- kstest(b,'alpha', (5/100)/Kmax);

    xlabel('CRBETA')
    ylabel('count')

    if rem(idx,25)==0
        drawnow
    end

    % pause(0.2)
end


passing_normal_crmsc = mean(normality_crmsc)
passing_normal_crbeta = mean(normality_crbeta)

%%
figure
subplot(2,1,1)

noise_freq_bins = 351:1:451;

hold on
for i = 1:numel(noise_freq_bins)
    fo = noise_freq_bins(i); 
    p1 = plot(cumsum(rolling_msc(fo,:)), 'Color', 0.8*[1 1 1],'Linewidth', 0.8);
end
for i = 1:numel(signal_freq_bins)
    fo = signal_freq_bins(i); 
    p2 = plot(cumsum(rolling_msc(fo,:)),'Linewidth', 1.2);
end

legend([p1,p2], 'f(ASSR)', 'f(Noise)', 'Location', 'northwest', 'Fontsize',14)
ylabel('Acumulated Magnitude-Squared Coherence [\mu V]')
grid on
% title(['MSC acumulada por janela, (',type,')',' mean(SNR) = ',num2str(noise_mean)])
hold off

%%
subplot(2,1,2)
if upperlim> size(rolling_msc,2)
    upperlim = size(rolling_msc,2);
end

K = Kmax;
aa = cumsum(rolling_msc(noise_freq_bins,:),2);
aa = aa(:,1:step:upperlim);
p1 = boxplot(aa,'Symbol','kx','OutlierSize',5,'Color', 0.7*[1 1 1]);
hold on
p2 = plot(quantile(aa,1-5/100), '-', 'Linewidth', 1.2, 'Color',[1 0.4 0.4]);

bb = cumsum(rolling_msc(signal_freq_bins,:),2);
bb = bb(:,1:step:upperlim);
p3 = boxplot(bb,'Symbol','bo','OutlierSize',5);
xlabel(' Time [100 ms]', 'Fontsize', 14)
ylabel(' Acumulated Magnitude-Squared Coherence [\mu V]')

%%
figure;
subplot(211)
stem(normality_crmsc, 'filled')
title(['Complement of decision from KS Test (',type,')'], 'Fontsize', 14)
ylabel('Is CR-MSC normal? [1- pval > 99% / K_{max}]', 'Fontsize', 14)
ylim([0,1.2])
grid on
subplot(212)
stem(normality_crbeta, 'filled')
ylabel('Is CR-\beta normal? [1- pval > 99% / K_{max}]', 'Fontsize', 14)
xlabel('Test index [K]','Fontsize', 14)
ylim([0,1.2])
grid on
end
end