% OBJ: gerar curva eficiencia de pareto para CGST MSC
% https://pt.wikipedia.org/wiki/Efici%C3%AAncia_de_Pareto

%useful viz plot:
% figure; stem(all_freqs,-2*log(Ps(all_freqs)))
%% Setup
clearvars; 
close all
clc

should_plot = 1;
%% Parametros
Kd = 88;           % numero total de testes a aplicar sequencialmente
FPd = 0.05;      % taxa de falso positivo desejado para o exame

Msample = 0; % best = 32 / 20 a 24
Mmin = 32; %floor(Msample*0.6);
Mstep = 32;
sum_step = Mstep;

NumT = 5e3; 
Nsinais = NumT;

int_inic =3;
max_int = 4;

FNRt = nan(11,max_int-int_inic+1, Kd);
TPRt = nan(11,max_int-int_inic+1, Kd);
TNRt = nan(11,max_int-int_inic+1, Kd);
FPRt = nan(11,max_int-int_inic+1, Kd);

for cont_int = int_inic:max_int
for cont_vol=1:11

% always plot last run
if cont_int ==max_int
    if cont_vol==11
        should_plot = 1;
    end
end



    cont_vol

% REAL: s frequências das portadoras para ambas os ouvidos foram as mesmas: 
% 500, 1000, 2000,4000 Hz, modulados, respectivamente, nas frequências 
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
fcInferior = 70; % 70
fcSuperior = Fs/2 -1; % 100
[b,a] = butter(8,[fcInferior/(Fs/2), fcSuperior/(Fs/2)]);
x = filter(b,a,x); 
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


TotalAlpha      = FPd;                        
Alpha_k         = ones(1,K)*(TotalAlpha/K);     
Gamma_k         = ((1-TotalAlpha)/K).*ones(1,K);
% Alpha_k = exp(10*linspace(TotalAlpha/K, TotalAlpha, K));
% Alpha_k = flip((TotalAlpha.*Alpha_k./sum(Alpha_k)));
% Alpha_k = (TotalAlpha.*Alpha_k./sum(Alpha_k));
% Gamma_k = exp(linspace((1-TotalAlpha)/K, 1-TotalAlpha, K));
% Gamma_k = flip((1-TotalAlpha).*Gamma_k./sum(Gamma_k));
% Gamma_k = (1-TotalAlpha).*Gamma_k./sum(Gamma_k);
% Gamma_k = 0.01.*ones(1,K);
% Gamma_k(end) = 1-sum(Gamma_k(1:end-1))-TotalAlpha;

% aThresholds  = zeros(1,K);
% gThresholds  = zeros(1,K);

% primeiro = SIGNALS(:,1:Msample);


% accp = zeros(NFFT/2+1, Msample);
% for k=Mmin+1:Msample
%     if k-M<=0
%         ind_inicial = 1;
%         ind_final = k;
%         recorte = SIGNALS(:,ind_inicial:ind_final);
%         Matual = k;
%     else
%         ind_final = k;
%         ind_inicial = ind_final-Mmin+1;
%         recorte = primeiro(:,ind_inicial:ind_final);
%         Matual = M;
%     end
% 
%     % ind_inicial = M*(k-1)+1;
%     % ind_final = ind_inicial+M-1;
%     accp(:,k) = msc_fft(recorte,Matual);
% end


% freq = 1:500;
% freq(signal_freq_bins) = [];
% h0 = sum(accp(noise_freq_bins,:),2);
% h0 = accp(freq,:);
% muk1=mean(h0);
% sigmak1=std(h0);

% primeiro = msc_fft(primeiro,Msample);
% muk1=mean(primeiro(noise_freq_bins));
% sigmak1=std(primeiro(noise_freq_bins));


Resolution      = (1/0.001);                  
Xvalues         = 0:1/Resolution:15;            
Null         	= betapdf(Xvalues, 1, M-1); %normpdf(Xvalues,muk1,sigmak1);% normpdf(Xvalues,0,0.05); %normpdf(Xvalues,0.05,0.3); 
Null            = Null/sum(Null);             	
Chi2_Norm       = Null/sum(Null);             	
k               = 1;                            
aThresholds(k)	= 1 - Alpha_k(k).^(1./(M-1));  % quantile(betarandn(1,1e5), 1-Alpha_k(1)); %
gThresholds(k)	= 1-(1- Gamma_k(k)).^(1./(M-1)); % quantile(randn(1,1e5), Gamma_k(1));
TruncInd_Ra      = round(aThresholds(k)*Resolution);                                         
TruncInd_Rg      = round(gThresholds(k)*Resolution);           

for k = 2:K
    % disp(k)
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


%% Plot thresholds
disp('')
disp('--------------------------------')
scr_siz = get(0,'ScreenSize') ;

if should_plot
    % floor([0.5*scr_siz(3)/2 0.6*scr_siz(4)/2 1.2*scr_siz(3)/2 1.2*scr_siz(4)/2])
    f=figure(1);
    f.Position =floor([100 100 1.5*scr_siz(3)/2 1.5*scr_siz(4)/2]);
    subplot(121)
    cv_max = max([aThresholds, gThresholds]);
    area([0,M],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha',0.1)
    hold on
    grid on
    area([M,K],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0 0.4470 0.7410], 'FaceAlpha',0.1)            
    
    dim = [0.14 0.75 0.1 0.1];
    str = {'Data collection', ['[M_{min}= ',num2str(M),']']};
    ta = annotation('textbox',dim,'String',str, ...
        'FitBoxToText','on', 'FontSize',12);
    ta.FaceAlpha = 0.2;
    ta.EdgeColor = [0.8500 0.3250 0.0980];  
    ta.BackgroundColor = [0.8500 0.3250 0.0980];  
    ta.Color = [.2 .2 .2]; 
    
    dim = [0.225 0.75 0.1 0.1];
    str = {'Test region'};
    tb = annotation('textbox',dim,'String',str, ...
        'FitBoxToText','on', 'FontSize',12);
    tb.FaceAlpha = 0.2;  
    tb.EdgeColor = [0 0.4470 0.7410]; 
    tb.BackgroundColor = [0 0.4470 0.7410]; 
    tb.Color = [.2 .2 .2]; 
end




%signal freqs
random_exam = zeros(numel(signal_freq_bins),K);
hit = zeros(1,K);
miss= zeros(1,K);

% c = 0.7*lines(numel(signal_freq_bins)); % sky, hot, turbo, parula, cool,
% spring
c = 0.9*gray(numel(signal_freq_bins));

for idx_freq = 1:numel(signal_freq_bins)
    
    flag = 1;
    for k=Mmin+1:K
        if k-M<=0
            ind_inicial = 1;
            ind_final = k*Mstep;
            recorte = SIGNALS(:,ind_inicial:ind_final);
            Matual = k;
        else
            ind_final = k;
            ind_inicial = ind_final-M+1;
            recorte = SIGNALS(:,ind_inicial:ind_final);
            Matual = M;
        end

        % ind_inicial = M*(k-1)+1;
        % ind_final = ind_inicial+M-1;
        Ps = msc_fft(recorte,Matual);
        random_exam(idx_freq,k) = Ps(signal_freq_bins(idx_freq));

        if sum(random_exam(idx_freq,Mmin:sum_step:k)) >= aThresholds(idx_freq) && flag && k>=Msample
            hit(k) = hit(k)+1;
            last_k = k;
            final_v = sum(random_exam(idx_freq,1:k));
            flag = 0;
            
            % break
        elseif sum(random_exam(idx_freq,Mmin:sum_step:k)) <= gThresholds(idx_freq) && flag && k>=Msample
            miss(k) = miss(k)+1;
            last_k = k;
            final_v = sum(random_exam(idx_freq,1:k));
            flag = 0;
            
            % break
        end
    end

    if should_plot
        
        plot(Mmin:sum_step:K, cumsum(random_exam(idx_freq,Mmin:sum_step:end)),'--', ...
            'LineWidth',1.5, 'Color', c(idx_freq,:))
        hold on 
        plot(Mmin:sum_step:K, cumsum(random_exam(idx_freq,Mmin:sum_step:end)),'.', ...
            'MarkerSize',12, 'Color', c(idx_freq,:))
    end
end

if should_plot
    subplot(121)
    colormap(c)
    cb = colorbar;
    cb.TickLabels = signal_freq_bins(round(linspace(1,numel(signal_freq_bins),11)));
    
    % p4 = plot(last_k, final_v,'o','MarskerSize',7, 'LineWidth',2,'Color',"#0072BD");
    
    p1 = plot(Mmin:K+Mmin-1, aThresholds, 'LineWidth',1.8, 'Color',"#77AC30"); %[0 0.4470 0.7410]);
    p2 = plot(Mmin:K+Mmin-1, gThresholds,'LineWidth',1.8, 'Color',"#A2142F");%[0.8500 0.3250 0.0980]);
    
    title(['Critical values for coherence-based early stopping exam with ',...
        num2str(round(100*FPd)),'% significance'], 'FontSize', 18)
    % legend([p1,p2],'Detection [ \alpha_k ]', 'Futility [ \gamma_k ]', ...
        % 'Fontsize', 15, 'Location', 'SouthEast', 'AutoUpdate', 'off')
    ylabel('\Sigma_1^k MSC  [summary statistic]', 'Fontsize', 20)
    xlabel('Exam k-th epoch [seconds]', 'Fontsize', 20)
    % c = 0.7*hot(numel(signal_freq_bins));
    
    xlim([0,K])
    ylim([0,1.05*cv_max])
    hold off
    drawnow

end


FNRs = miss/numel(signal_freq_bins);
TPRs = hit/numel(signal_freq_bins);
FNRt(cont_vol,cont_int-int_inic+1,1:K) = FNRs;
TPRt(cont_vol,cont_int-int_inic+1,1:K) = TPRs;

if should_plot
%noise freqs
c2 = 0.9*gray(numel(noise_freq_bins));
subplot(122)
cv_max = max([aThresholds, gThresholds]);
    area([0,M],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha',0.1)
    hold on
    grid on
    area([M,K],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0 0.4470 0.7410], 'FaceAlpha',0.1)            
    
    dim = [0.14 0.75 0.1 0.1];
    str = {'Data collection', ['[M_{min}= ',num2str(M),']']};
    ta = annotation('textbox',dim,'String',str, ...
        'FitBoxToText','on', 'FontSize',12);
    ta.FaceAlpha = 0.2;
    ta.EdgeColor = [0.8500 0.3250 0.0980];  
    ta.BackgroundColor = [0.8500 0.3250 0.0980];  
    ta.Color = [.2 .2 .2]; 
    
    dim = [0.225 0.75 0.1 0.1];
    str = {'Test region'};
    tb = annotation('textbox',dim,'String',str, ...
        'FitBoxToText','on', 'FontSize',12);
    tb.FaceAlpha = 0.2;  
    tb.EdgeColor = [0 0.4470 0.7410]; 
    tb.BackgroundColor = [0 0.4470 0.7410]; 
    tb.Color = [.2 .2 .2]; 
end

random_exam = zeros(numel(noise_freq_bins),K);
hit = zeros(1,K);
miss= zeros(1,K);

for idx_freq = 1:numel(noise_freq_bins)
    flag = 1;
    for k=Mmin+1:K
        if k-M<=0
            ind_inicial = 1;
            ind_final = k;
            recorte = SIGNALS(:,ind_inicial:ind_final);
            Matual = k;
        else
            ind_final = k;
            ind_inicial = ind_final-M+1;
            recorte = SIGNALS(:,ind_inicial:ind_final);
            Matual = M;
        end
    
        % ind_inicial = M*(k-1)+1;
        % ind_final = ind_inicial+M-1;
        Ps = msc_fft(recorte,Matual);
        random_exam(idx_freq,k) = Ps(noise_freq_bins(idx_freq));
    
        if sum(random_exam(idx_freq,Mmin:k)) >= aThresholds(k) && flag && k>=Msample
            hit(k) = hit(k)+1;
            last_k = k;
            final_v = sum(random_exam(idx_freq,1:sum_step:k));
            flag = 0;
            
            % break
        elseif sum(random_exam(idx_freq,Mmin:k)) <= gThresholds(k) && flag && k>=Msample
            miss(k) = miss(k)+1;
            last_k = k;
            final_v = sum(random_exam(idx_freq,1:sum_step:k));
            flag = 0;
            
            % break
        end
    
        
    end
    if should_plot
            plot(Mmin:sum_step:K, cumsum(random_exam(idx_freq,Mmin:sum_step:end)),'--', ...
                'LineWidth',1.5, 'Color', c2(idx_freq,:))
            hold on
            plot(Mmin:sum_step:K, cumsum(random_exam(idx_freq,Mmin:sum_step:end)),'.', ...
                'MarkerSize',12, 'Color', c2(idx_freq,:))
        end
end

if should_plot
    colormap(c2)
    cb = colorbar;
    cb.TickLabels = noise_freq_bins(round(linspace(1,numel(noise_freq_bins),11)));
    
    % p4 = plot(last_k, final_v,'o','MarskerSize',7, 'LineWidth',2,'Color',"#0072BD");
    
    p1 = plot(Mmin:K+Mmin-1, aThresholds, 'LineWidth',1.8, 'Color',"#77AC30"); %[0 0.4470 0.7410]);
    hold on 
    p2 = plot(Mmin:K+Mmin-1, gThresholds,'LineWidth',1.8, 'Color',"#A2142F");%[0.8500 0.3250 0.0980]);
    
    % title(['Critical values for coherence-based early stopping exam with ',...
    %     num2str(round(100*FPd)),'% significance'], 'FontSize', 18)
    legend([p1,p2],'Detection [ \alpha_k ]', 'Futility [ \gamma_k ]', ...
        'Fontsize', 15, 'Location', 'NorthWest', 'AutoUpdate', 'off')
    ylabel('\Sigma_1^k MSC  [summary statistic]', 'Fontsize', 20)
    xlabel('Exam k-th epoch [seconds]', 'Fontsize', 20)
    % c = 0.7*hot(numel(signal_freq_bins));
    
    xlim([0,K])
    ylim([0,1.05*cv_max])
    hold off

    drawnow
    pause(0.1)
end



TNRs = miss/numel(noise_freq_bins);
FPRs = hit/numel(noise_freq_bins);
TNRt(cont_vol,cont_int-int_inic+1,1:K) = TNRs;
FPRt(cont_vol,cont_int-int_inic+1,1:K) = FPRs;

end

mTP = mean(TPRt,3,'omitnan')
mFP = mean(FPRt,3,'omitnan')
mFN = mean(FNRt,3,'omitnan')
mTN = mean(TNRt,3,'omitnan')


end
%%
% Taxa de deteccao media por estágio para cada intensidade:
pd_k =  squeeze(mean(TPRt,1, 'omitnan')); % intensidade = linhas, estagio = colunas

% Taxa de deteccao total por intensidade:
PD = round(100*sum(pd_k,2),2);

% similarmente para outros:
fp_k = squeeze(mean(FPRt,1,'omitnan'));
tn_k = squeeze(mean(TNRt,1,'omitnan'));
fn_k = squeeze(mean(FNRt,1,'omitnan'));

FP = round(100*sum(fp_k,2),2);
TN = round(100*sum(tn_k,2),2);
FN = round(100*sum(fn_k,2),2);

% mostrar resultados:
Resultados = table(cell2mat(Vintensidade(int_inic:max_int)), PD,FP, TN, FN, ...
    'VariableNames',{'Intensidade [dB]', ...
    'PD [%]', 'FP \alpha [%]', ...
    'TN [%]', 'FN [%] \beta'});
fprintf('\n\n\n')
disp(Resultados)


% fp_final = mean(FPRt,'all')
% pd_final = mean(TPRt,'all')
% mean(FNRt,'all')
% mean(TNRt,'all')

%% mostrar resultados 
hit = 100*pd_k;
miss = 100*fp_k;

%clc
%1 - Falsos Positivo  
figure 
plot(hit,'.k','MarkerSize',10)
hold on 
plot([0 size(hit,2)],[hit(end) hit(end)], ':r','LineWidth',2)
ylabel('Taxa de Detecção','fontsize',12)
grid on


%2 - Taxa de detecção 
FP_desejado = sum(Alpha_k);
figure 
plot(miss,'.k','MarkerSize',10)
hold on 
plot([0 size(miss,2)],[100*FP_desejado 100*FPd], ':r','LineWidth',2)
ylabel('Falso Positivo','fontsize',12)

figure
boxplot(miss)
hold on
plot([0 size(miss,2)],[100*FP_desejado 100*FPd], ':r','LineWidth',2)
