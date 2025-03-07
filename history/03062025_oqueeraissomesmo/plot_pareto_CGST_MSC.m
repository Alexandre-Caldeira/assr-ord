% OBJ: gerar curva eficiencia de pareto para CGST MSC
% https://pt.wikipedia.org/wiki/Efici%C3%AAncia_de_Pareto

%% Setup
clearvars; close all; clc

%% Parametros
K = 13;           % numero total de testes a aplicar sequencialmente
FPd = 0.05;      % taxa de falso positivo desejado para o exame

NumT = 5e3;
FS = 1000;  % frequencia de amostragem 
SFREQ = 88; % frequencia de estimulacao 
NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
Nsinais = NumT;

% snr_vec = 0:-1:-30;
snr_vec = 0:-1:-30;

Stage_FPRs  = zeros(numel(snr_vec),K);        % stage-wise TPRs
Stage_TNRs  = zeros(numel(snr_vec),K);        % stage-wise FNRs
FPR         = zeros(1,numel(snr_vec));        % total TPR
TNR         = zeros(1,numel(snr_vec));        % total FNR

Stage_TPRs  = zeros(numel(snr_vec),K);        % stage-wise TPRs
Stage_FNRs  =  zeros(numel(snr_vec),K);       % stage-wise FNRs
TPR         = zeros(1,numel(snr_vec));        % total TPR
FNR         = zeros(1,numel(snr_vec));        % total FNR


%%%%%%%%%%%%
cont_vol = 2
cont_int = 5
signal_freq_bins =  [82   90    84    86    88    90    92    94    96];
noise_freq_bins = [300 400 500] %round(signal_freq_bins.*exp(1)/2)+5;
%%%%%%%%%%%%
%% Carrega dados

% SET THIS PATH:
path = 'C:\PPGEE\SBEB_CBA_24\ASSR - Coleta OFFLINE';
addpath(path)
%vetor dos voluntários 
Vvoluntario = {'Abdon';'Ana';'BBB';'Colatina';'Erick';'Luciana';...
    'Sombra';'Quenaz';'Vinicius';'Sacola';'Wreikson'}; 

%vetor da intensidade 
Vintensidade = {'70';'60';'50';'32';'30'}; 
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

fcInferior = 70;
fcSuperior = 110;

d2 = designfilt('bandpassiir','FilterOrder',2, ...
    'HalfPowerFrequency1',fcInferior,'HalfPowerFrequency2',fcSuperior, ...
    'DesignMethod','butter','SampleRate',Fs);

x = filtfilt(d2,x);

% x = filter(b,a,x); 
% excluir os dois primeiros segundos do inicio da coleta 
x(:,1:2,:) =[]; 

%encontrar o valor máximo por canal 
Vmax = squeeze(max(max(abs(x)),[],3));
ind = Vmax>remoc;

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
pos_eletrodo= 1;
xmedia = x(:,~ind,pos_eletrodo);
% xmedia = x(:,:,pos_eletrodo);

freq=0:(Fs-1); 

% xmedia = xmedia./std(xmedia);

SIGNALS = fft(xmedia(:,:))*2/nfft*1e9;
FS = Fs;
NFFT = nfft;
SIGNALS = SIGNALS(1:floor(end/2)+1,:); % only half the FFT spectrum is valid
f = FS/2*linspace(0,1,NFFT/2+1)'; % only half the FFT spectrum is valid
max_length = size(SIGNALS,2);

all_freqs = [signal_freq_bins noise_freq_bins];
% 
% for idx_epoch = 1:max_length 
% 
%     % para cada frequencia em f
%     for idx_f = 1:1:numel(f)
% 
% 
%         M = idx_epoch;
%         if M>32
%             M=32;
% 
%             X_atual = SIGNALS(:,idx_epoch-M+1:idx_epoch);
% 
%         else
% 
%             X_atual = SIGNALS(:,1:idx_epoch);             
% 
%         end
% 
%         % GFT(idx_f,idx_epoch) = sum(abs(X_atual(idx_f,1:M,:)).^2)./...
%         %             (sum(abs(X_atual(idx_f,1:M,:)).^2)+sum(abs(X_atual(numel(f),1:M,:)).^2));
%         % 
%         % c1_csm = (sum(cos(angle(X_atual(idx_f,:))))./M).^2;
%         % c2_csm = (sum(sin(angle(X_atual(idx_f,:))))./M).^2;
%         % CSM(idx_f,idx_epoch) = c1_csm+c2_csm;
% 
%         num_msc = abs(sum(X_atual(idx_f,:)))^2;
%         den_msc= M*sum(abs(X_atual(idx_f,:)).^2);
%         MSC(idx_f,idx_epoch) = num_msc/den_msc;
%     end
% end
% 
% % c_states(:,3,:) = log10(MSC(all_freqs,:));
% MSC = MSC(all_freqs,2:end);

%% Calc thresholds
% K               = 5;                            
M = floor(size(SIGNALS,2)/K);
TotalAlpha      = FPd;                        
% Alpha_k         = ones(1,K)*(TotalAlpha/K);     
% Alpha_k = linspace(TotalAlpha/2*K, TotalAlpha, K);
% Alpha_k = TotalAlpha.*Alpha_k./sum(Alpha_k);
% % Gamma_k         = ((1-TotalAlpha)/K).*ones(1,K);
% Gamma_k = linspace((1-TotalAlpha)/2*K, 1-TotalAlpha, K);
% Gamma_k = (1-TotalAlpha).*Gamma_k./sum(Gamma_k);
% Gamma_k = [0.0001 0.0001 0.01 0.9398];

% Alpha_k         = ones(1,K)*(TotalAlpha/K);     
Alpha_k = exp(10*linspace(TotalAlpha/K, TotalAlpha, K));
Alpha_k = flip((TotalAlpha.*Alpha_k./sum(Alpha_k)));
% Gamma_k         = ((1-TotalAlpha)/K).*ones(1,K);
Gamma_k = exp(linspace((1-TotalAlpha)/K, 1-TotalAlpha, K));
Gamma_k = flip((1-TotalAlpha).*Gamma_k./sum(Gamma_k));

Resolution      = (1/0.0001);                  
Xvalues         = 0:1/Resolution:1;            
Null         	= betapdf(Xvalues, 1, M-1); 
Null            = Null/sum(Null);             	
Chi2_Norm       = Null/sum(Null);             	
k               = 1;                            
aThresholds(k)	= 1 - Alpha_k(k).^(1./(M-1));  
gThresholds(k)	= 1-(1- Gamma_k(k)).^(1./(M-1));
TruncInd_Ra      = round(aThresholds(k)*Resolution);                                         
TruncInd_Rg      = round(gThresholds(k)*Resolution);           

for k = 2:K
    NullTrunc                   = Null;                                                     
    NullTrunc(TruncInd_Ra:end)  = zeros(1, length(NullTrunc(TruncInd_Ra:end)));              
    NullTrunc(1:TruncInd_Rg)    = zeros(1, length(NullTrunc(1:TruncInd_Rg)));
    Null2                       = conv(Chi2_Norm, NullTrunc);                              
    Null2                       = Null2 / (sum(Null2) / (1 - sum(Gamma_k(1:(k-1))) - sum(Alpha_k(1:(k-1)))));
    TruncInd_Ra                 = findIndex(Null2, sum(Null2) - Alpha_k(k));            
    aThresholds(k)              = TruncInd_Ra/Resolution;                                    
    TruncInd_Rg                 = findIndex(Null2, Gamma_k(k), 1);
    gThresholds(k)              = TruncInd_Rg/Resolution;
    Null                        = Null2;
end

%%
SFREQ = 1;

count_tests = zeros(1,K);
Stage_FPRs  = zeros(1,K);        % stage-wise TPRs
Stage_TNRs  = zeros(1,K);        % stage-wise FNRs
FPR         = zeros(1,1);        % total TPR
TNR         = zeros(1,1);        % total FNR

Stage_TPRs  = zeros(1,K);        % stage-wise TPRs
Stage_FNRs  =  zeros(1,K);       % stage-wise FNRs
TPR         = zeros(1,1);        % total TPR
FNR         = zeros(1,1);        % total FNR

% get TP, FN
count = 0;
TP      = zeros(1, K);      % number of false-positives
FN      = zeros(1, K);      % number of true-negatives

NumJ = size(SIGNALS,2); % numero de janelas total (1 teste por janela)
% for ti=1:NumJ    
count = count+1;
Plog = zeros(1,K);

for k=1:K        
    count_tests(k) = count_tests(k)+1;
    ind_inicial = M*(k-1)+1;
    ind_final = ind_inicial+M-1;

    Ps = msc_fft(SIGNALS(:,ind_inicial:ind_final),M);
    % Ps = MSC(:,ti);
    Plog(k) = Ps(SFREQ);

    if sum(Plog(1:k)) >= aThresholds(k)
        TP(k) = TP(k)+1;
        break

    elseif sum(Plog(1:k)) <= gThresholds(k)
        FN(k) = FN(k)+1;
        break

    end

end
% end

% get FP, TN
count = 0;
FP      = zeros(1, K);      % number of false-positives
TN      = zeros(1, K);      % number of true-negatives

% for ti=1:NumJ
    % Ps          = msc_fft(S5(:,:,ti),M);
    % Plog        = Ps;           
count = count+1;
% check for rejections
Plog = zeros(1,K);




for k=1:K        

    ind_inicial = M*(k-1)+1;
    ind_final = ind_inicial+M-1;
    Ps = msc_fft(SIGNALS(all_freqs,ind_inicial:ind_final),M);
    Plog(k) = Ps(end);
    % Ps = msc_fft(S5(:,ind_inicial:ind_final,ti),M);
    % Ps = MSC(:,ti);
    % Plog(k) = Ps(end);
    

    if sum(Plog(1:k)) >= aThresholds(k)
        FP(k) = FP(k)+1;
        break

    elseif sum(Plog(1:k)) <= gThresholds(k)
        TN(k) = TN(k)+1;
        break

    end

end
% end

Stage_FPRs   = FP/count               % stage-wise TPRs
Stage_TNRs   = TN/count               % stage-wise FNRs
Stage_TPRs  = TP/count               % stage-wise TPRs
Stage_FNRs  = FN/count               % stage-wise FNRs
TPR         = sum(TP) / count        % total TPR
FNR        = sum(FN) / count        % total FNR
FPR         = sum(FP) / count        % total TPR
TNR         = sum(TN) / count        % total FNR


%% Plot thresholds
disp('')
disp('--------------------------------')

figure %(1)

cv_max = max([aThresholds, gThresholds]);
area([0,M],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha',0.3)
hold on
area([M,K*M],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0 0.4470 0.7410], 'FaceAlpha',0.3)            

dim = [0.15 0.55 0.15 0.3];
str = {'Data collection', ['[M_{min}= ',num2str(M),']']};
ta = annotation('textbox',dim,'String',str, ...
    'FitBoxToText','on', 'FontSize',12);
ta.FaceAlpha = 0.3;
ta.EdgeColor = [0.8500 0.3250 0.0980];  
ta.BackgroundColor = [0.8500 0.3250 0.0980];  
ta.Color = [.2 .2 .2]; 

dim = [0.3 0.55 0.3 0.3];
str = {'Test region'};
tb = annotation('textbox',dim,'String',str, ...
    'FitBoxToText','on', 'FontSize',12);
tb.FaceAlpha = 0.3;  
tb.EdgeColor = [0 0.4470 0.7410]; 
tb.BackgroundColor = [0 0.4470 0.7410]; 
tb.Color = [.2 .2 .2]; 

% random_exam = cumsum(betarnd(1,M-1,1,K));
random_exam = zeros(1,K);
flag = 1;
for k=1:K
    ind_inicial = M*(k-1)+1;
    ind_final = ind_inicial+M-1;
    Ps = msc_fft(SIGNALS(:,ind_inicial:ind_final),M);
    random_exam(k) = Ps(signal_freq_bins(randi(numel(signal_freq_bins),1)));
    if sum(random_exam(1:k)) >= aThresholds(k) && flag
        % random_exam(k+1:K) = NaN;
        last_k = k;
        final_v = sum(random_exam(1:k));
        flag = 0;
        % break
    elseif sum(random_exam(1:k)) <= gThresholds(k) && flag
        % random_exam(k+1:K) = NaN;
        last_k = k;
        final_v = sum(random_exam(1:k));
        flag = 0;
        % break
    end
end

p1 = plot(M:M:M*K, aThresholds, 'LineWidth',1.2, 'Color',"#77AC30"); %[0 0.4470 0.7410]);
p2 = plot(M:M:M*K, gThresholds,'LineWidth',1.2, 'Color',"#A2142F");%[0.8500 0.3250 0.0980]);
p3 = plot(M:M:M*K, cumsum(random_exam),'LineWidth',1.2, 'Color',"#0072BD");
p4 = plot(last_k, final_v,'o','MarkerSize',7, 'LineWidth',2,'Color',"#0072BD");
grid on

title('Critical values for coherence-based early stopping exam with 1% significance', 'FontSize', 15)
legend([p1,p2,p3,p4],'\alpha_k (Detection)', '\gamma_k (Futility)','Exam','Early Stop', ...
    'Fontsize', 15, 'Location', 'SouthEast')
ylabel('\Sigma_1^k MSC  [summary statistic]', 'Fontsize', 20)
xlabel('Exam k-th epoch [seconds]', 'Fontsize', 20)

hold off

xlim([0,K*M])
ylim([0,1.05*cv_max])
drawnow













