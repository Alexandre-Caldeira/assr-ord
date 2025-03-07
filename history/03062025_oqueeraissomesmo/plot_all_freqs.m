% OBJ: gerar curva eficiencia de pareto para CGST MSC
% https://pt.wikipedia.org/wiki/Efici%C3%AAncia_de_Pareto

%% Setup
clearvars; 
% close all;
clc

%% Parametros
Kd = 10;           % numero total de testes a aplicar sequencialmente
FPd = 0.05;      % taxa de falso positivo desejado para o exame
Mmin = 32;

NumT = 5e3; 
Nsinais = NumT;

% snr_vec = 0:-1:-30;
snr_vec = 0:-1:-30;
for cont_vol=1:11
%%%%%%%%%%%%
% cont_vol = 1;
% cont_int = 5;

% cont_vol = 11

% cont_int = 3;
for cont_int = 2:5
signal_freq_bins =  [82  84  86    88    90    92    94    96];
noise_freq_bins = 400:5:441;%300+round(signal_freq_bins.*exp(1)/2)+5; %[300 400 500]; %
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
fcInferior = 70;
fcSuperior = 100;


d1 = designfilt('bandstopiir','FilterOrder',4, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);

d2 = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',fcInferior,'HalfPowerFrequency2',fcSuperior, ...
    'DesignMethod','butter','SampleRate',Fs);

x = filtfilt(d1,x);
x = filtfilt(d2,x);

%colatina
fcInferior = 70;
fcSuperior = 100;
[b,a] = butter(2,[fcInferior/(Fs/2), fcSuperior/(Fs/2)]);
x = filter(b,a,x); 
% excluir os dois primeiros segundos do inicio da coleta 
x(:,1:2,:) =[]; 

% %encontrar o valor máximo por canal 
% Vmax = squeeze(max(max(abs(x)),[],3));
% ind = Vmax>remoc;
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
pos_eletrodo= 1;
% xmedia = x(:,~ind,pos_eletrodo);
xmedia = x(:,:,pos_eletrodo);

freq=0:(Fs-1); 

% % xmedia = xmedia./std(xmedia);

SIGNALS = fft(xmedia(:,:))*2/nfft*1e9;
FS = size(SIGNALS,1);
NFFT = size(SIGNALS,1);
SIGNALS = SIGNALS(1:floor(end/2)+1,:); % only half the FFT spectrum is valid
f = FS/2*linspace(0,1,NFFT/2+1)'; % only half the FFT spectrum is valid
max_length = size(SIGNALS,2);

all_freqs = [signal_freq_bins noise_freq_bins];
SIGNALS = SIGNALS(:,3:end);
SIGNALS = SIGNALS./std(SIGNALS,0,2);
% SIGNALS = SIGNALS- mean(SIGNALS,2);

%% Calc thresholds
                     
M = floor(size(SIGNALS,2)/Kd);
if M<Mmin
    if floor(size(SIGNALS,2)/Mmin)<=1
        % exame_tem_apenasfloor(size(SIGNALS,2)/Kd))
        
        % error
        disp(['Exame curto demais para teste sequencial com MSC! ', ...
           num2str(size(SIGNALS,2)),' janelas totais.' ])
        break
    end  

    K = floor(size(SIGNALS,2)/Mmin)
    M=Mmin
else
    K = Kd
    M
end


TotalAlpha      = FPd;                        
Alpha_k         = ones(1,K)*(TotalAlpha/K);     
% Alpha_k = linspace(TotalAlpha/2*K, TotalAlpha, K);
% Alpha_k = TotalAlpha.*Alpha_k./sum(Alpha_k);
Gamma_k         = ((1-TotalAlpha)/K).*ones(1,K);
% Gamma_k = linspace((1-TotalAlpha)/2*K, 1-TotalAlpha, K);
% Gamma_k = (1-TotalAlpha).*Gamma_k./sum(Gamma_k);
% Gamma_k = [0.0001 0.0001 0.01 0.9398];
% Gamma_k = 0.01.*ones(1,K);
% Gamma_k(end) = 1-sum(Gamma_k(1:end-1))-TotalAlpha;

% Alpha_k         = ones(1,K)*(TotalAlpha/K);     
% Alpha_k = exp(100*linspace(TotalAlpha/K, TotalAlpha, K));
% Alpha_k = flip((TotalAlpha.*Alpha_k./sum(Alpha_k)));
% % Gamma_k         = ((1-TotalAlpha)/K).*ones(1,K);
% Gamma_k = exp(linspace((1-TotalAlpha)/K, 1-TotalAlpha, K));
% Gamma_k = flip((1-TotalAlpha).*Gamma_k./sum(Gamma_k));

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

% Stage_FPRs   = FP/count               % stage-wise TPRs
% Stage_TNRs   = TN/count               % stage-wise FNRs
% Stage_TPRs  = TP/count               % stage-wise TPRs
% Stage_FNRs  = FN/count               % stage-wise FNRs
% TPR         = sum(TP) / count        % total TPR
% FNR        = sum(FN) / count        % total FNR
% FPR         = sum(FP) / count        % total TPR
% TNR         = sum(TN) / count        % total FNR


%% Plot thresholds
disp('')
disp('--------------------------------')


figure(1)

cv_max = max([aThresholds, gThresholds]);
area([0,M],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha',0.1)
hold on
grid on
area([M,K*M],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0 0.4470 0.7410], 'FaceAlpha',0.1)            

dim = [0.15 0.55 0.15 0.3];
str = {'Data collection', ['[M_{min}= ',num2str(M),']']};
ta = annotation('textbox',dim,'String',str, ...
    'FitBoxToText','on', 'FontSize',12);
ta.FaceAlpha = 0.2;
ta.EdgeColor = [0.8500 0.3250 0.0980];  
ta.BackgroundColor = [0.8500 0.3250 0.0980];  
ta.Color = [.2 .2 .2]; 

dim = [0.3 0.55 0.3 0.3];
str = {'Test region'};
tb = annotation('textbox',dim,'String',str, ...
    'FitBoxToText','on', 'FontSize',12);
tb.FaceAlpha = 0.2;  
tb.EdgeColor = [0 0.4470 0.7410]; 
tb.BackgroundColor = [0 0.4470 0.7410]; 
tb.Color = [.2 .2 .2]; 


%signal freqs
random_exam = zeros(numel(signal_freq_bins),K);
hit = zeros(1,K);
miss= zeros(1,K);

% c = 0.7*lines(numel(signal_freq_bins)); % sky, hot, turbo, parula, cool,
% spring
c = 0.7*gray(numel(signal_freq_bins));

for idx_freq = 1:numel(signal_freq_bins)
    
    flag = 1;
    for k=1:K
        ind_inicial = M*(k-1)+1;
        ind_final = ind_inicial+M-1;
        Ps = msc_fft(SIGNALS(:,ind_inicial:ind_final),M);
        random_exam(idx_freq,k) = Ps(signal_freq_bins(idx_freq));

        if sum(random_exam(idx_freq,1:k)) >= aThresholds(k) && flag
            hit(k) = hit(k)+1;
            last_k = k;
            final_v = sum(random_exam(idx_freq,1:k));
            flag = 0;
            
            % break
        elseif sum(random_exam(idx_freq,1:k)) <= gThresholds(k) && flag
            miss(k) = miss(k)+1;
            last_k = k;
            final_v = sum(random_exam(idx_freq,1:k));
            flag = 0;
            
            % break
        end
    end
    plot(M:M:M*K, cumsum(random_exam(idx_freq,:)),'--', ...
        'LineWidth',0.8, 'Color', c(idx_freq,:))
    plot(M:M:M*K, cumsum(random_exam(idx_freq,:)),'.', ...
        'MarkerSize',3, 'Color', c(idx_freq,:))
end

colormap(c)
cb = colorbar;
cb.TickLabels = signal_freq_bins(round(linspace(1,numel(signal_freq_bins),11)));

% p4 = plot(last_k, final_v,'o','MarskerSize',7, 'LineWidth',2,'Color',"#0072BD");

p1 = plot(M:M:M*K, aThresholds, 'LineWidth',1.2, 'Color',"#77AC30"); %[0 0.4470 0.7410]);
p2 = plot(M:M:M*K, gThresholds,'LineWidth',1.2, 'Color',"#A2142F");%[0.8500 0.3250 0.0980]);

title(['Critical values for coherence-based early stopping exam with ',...
    num2str(round(100*FPd)),'% significance'], 'FontSize', 15)
legend([p1,p2],'\alpha_k (Detection)', '\gamma_k (Futility)', ...
    'Fontsize', 15, 'Location', 'SouthEast', 'AutoUpdate', 'off')
ylabel('\Sigma_1^k MSC  [summary statistic]', 'Fontsize', 20)
xlabel('Exam k-th epoch [seconds]', 'Fontsize', 20)
% c = 0.7*hot(numel(signal_freq_bins));

xlim([0,K*M])
ylim([0,1.05*cv_max])
hold off
% drawnow


FNRs = miss/numel(signal_freq_bins);
TPRs = hit/numel(signal_freq_bins);
FNRt(cont_vol,cont_int) = sum(FNRs);
TPRt(cont_vol,cont_int) = sum(TPRs);


%noise freqs
random_exam = zeros(numel(signal_freq_bins),K);
hit = zeros(1,K);
miss= zeros(1,K);

for idx_freq = 1:numel(noise_freq_bins)
    
    flag = 1;
    for k=1:K
        ind_inicial = M*(k-1)+1;
        ind_final = ind_inicial+M-1;
        Ps = msc_fft(SIGNALS(:,ind_inicial:ind_final),M);
        random_exam(idx_freq,k) = Ps(noise_freq_bins(idx_freq));

        if sum(random_exam(idx_freq,1:k)) >= aThresholds(k) && flag
            hit(k) = hit(k)+1;
            last_k = k;
            final_v = sum(random_exam(idx_freq,1:k));
            flag = 0;
            
            % break
        elseif sum(random_exam(idx_freq,1:k)) <= gThresholds(k) && flag
            miss(k) = miss(k)+1;
            last_k = k;
            final_v = sum(random_exam(idx_freq,1:k));
            flag = 0;
            
            % break
        end
    end
end
TNRs = miss/numel(noise_freq_bins);
FPRs = hit/numel(noise_freq_bins);
TNRt(cont_vol,cont_int) = sum(TNRs);
FPRt(cont_vol,cont_int) = sum(FPRs);


end
end
%%
mTP = mean(TPRt)
mFP = mean(FPRt)
mFN = mean(FNRt)
mTN = mean(TNRt)




