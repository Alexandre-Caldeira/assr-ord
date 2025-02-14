% OBJ: gerar curva eficiencia de pareto para CGST MSC
% https://pt.wikipedia.org/wiki/Efici%C3%AAncia_de_Pareto

%useful viz plot:
% figure; stem(all_freqs,-2*log(Ps(all_freqs)))
%% Setup
clearvars; 
% close all;
clc

%% Parametros
Kd = 8;           % numero total de testes a aplicar sequencialmente
FPd = 0.05;      % taxa de falso positivo desejado para o exame
Mmin = 32;

NumT = 5e3; 
Nsinais = NumT;

for cont_vol=1:11
    
    max_int = 5;
    int_inic = 3;
    
    for cont_int = int_inic:max_int
        
        % REAL: s frequências das portadoras para ambas os ouvidos foram as mesmas: 
        % 500, 1000, 2000,4000 Hz, modulados, respectivamente, nas frequências 
        % 81, 85, 89 e 93 Hz, para o ouvido direito, e 83, 87, 91 e 95 Hz, 
        % para o ouvido esquerdo . - pg 52 tese colatina
        
        signal_freq_bins =  [82 84 86 88 90 92 94 96];
        noise_freq_bins = 351:1:451;

        %%%%%%%%%%%%
        %% Carrega dados
        
        % Set path to data:
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
        
        nfft = Fs; %1segundo de sinal 
         
        %retirar componente DC por janela (fiz isso pq no processamento em
        %tempo real é por janela)
        x = x - repmat(mean(x),nfft,1); %tirar a média de cada de cada trecho - devido a remoção
        

        % Em seguida, os sinais foram filtrados digitalmente por um filtro
        % Butterworth passa-banda de oitava ordem com frequência de corte de 1 Hz 
        % acima e abaixo da frequência de modulação (Antunese Felix, 2019).
        % - pg 53, tese colatina
        fcInferior = min(signal_freq_bins)-1;
        fcSuperior = max(signal_freq_bins)+1;

        [b,a] = butter(8,[fcInferior/(Fs/2), fcSuperior/(Fs/2)]);
        x = filter(b,a,x); 
        % excluir os dois primeiros segundos do inicio da coleta 
        x(:,1:2,:) =[]; 
        
        % Encontrar o valor máximo por canal 
        Vmax = squeeze(max(max(abs(x)),[],3));
      
        pos_eletrodo    = 1;
        xmedia = x(:,~ind,pos_eletrodo);
        
        SIGNALS = fft(xmedia,Fs); %*2/nfft*1e9;
        FS = Fs;
        NFFT = Fs;
        SIGNALS = SIGNALS(1:floor(end/2)+1,:); % only half the FFT spectrum is valid
        
        all_freqs = [signal_freq_bins noise_freq_bins];

        SIGNALS = normalize(SIGNALS);      
        
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
        Gamma_k         = ((1-TotalAlpha)/K).*ones(1,K);
        
        aThresholds  = zeros(1,K);
        gThresholds  = zeros(1,K);
        
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
        
        
        %% Plot thresholds
        disp('')
        disp('--------------------------------')
        should_plot = 1;
        
        if should_plot
            figure(1)
            
            cv_max = max([aThresholds, gThresholds]);
            area([0,M],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0.8500 0.3250 0.0980], 'FaceAlpha',0.1)
            hold on
            grid on
            area([M,K*M],[1.2*cv_max,1.2*cv_max], 'FaceColor',[0 0.4470 0.7410], 'FaceAlpha',0.1)            
            
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
        
        c = 0.9*gray(numel(signal_freq_bins));
        
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
                    
                    % continue
                elseif sum(random_exam(idx_freq,1:k)) <= gThresholds(k) && flag
                    miss(k) = miss(k)+1;
                    last_k = k;
                    final_v = sum(random_exam(idx_freq,1:k));
                    flag = 0;
                    
                    % continue
                end
            end
        
            if should_plot
                plot(M:M:M*K, cumsum(random_exam(idx_freq,:)),'--', ...
                    'LineWidth',1.5, 'Color', c(idx_freq,:))
                plot(M:M:M*K, cumsum(random_exam(idx_freq,:)),'.', ...
                    'MarkerSize',12, 'Color', c(idx_freq,:))
            end
        end
        
        if should_plot
            colormap(c)
            cb = colorbar;
            cb.TickLabels = signal_freq_bins(round(linspace(1,numel(signal_freq_bins),11)));

            p1 = plot(M:M:M*K, aThresholds, 'LineWidth',1.8, 'Color',"#77AC30"); %[0 0.4470 0.7410]);
            p2 = plot(M:M:M*K, gThresholds,'LineWidth',1.8, 'Color',"#A2142F");%[0.8500 0.3250 0.0980]);
            
            title(['Critical values for coherence-based early stopping exam with ',...
                num2str(round(100*FPd)),'% significance'], 'FontSize', 18)
            legend([p1,p2],'Detection [ \alpha_k ]', 'Futility [ \gamma_k ]', ...
                'Fontsize', 15, 'Location', 'SouthEast', 'AutoUpdate', 'off')
            ylabel('\Sigma_1^k MSC  [summary statistic]', 'Fontsize', 20)
            xlabel('Exam k-th epoch [seconds]', 'Fontsize', 20)

            xlim([0,K*M])
            ylim([0,1.05*cv_max])
            hold off
        end

        
        FNRs = miss/numel(signal_freq_bins);
        TPRs = hit/numel(signal_freq_bins);
        FNRt(cont_vol,cont_int-int_inic+1,:) = FNRs;
        TPRt(cont_vol,cont_int-int_inic+1,:) = TPRs;
            
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
        TNRt(cont_vol,cont_int-int_inic+1,:) = TNRs;
        FPRt(cont_vol,cont_int-int_inic+1,:) = FPRs;
        
    end
    
    % mTP = mean(TPRt)
    % mFP = mean(FPRt)
    % mFN = mean(FNRt)
    % mTN = mean(TNRt)


end
%%
% Taxa de deteccao media por estágio para cada intensidade:
pd_k =  squeeze(mean(TPRt,1)); % intensidade = linhas, estagio = colunas

% Taxa de deteccao total por intensidade:
PD = round(100*sum(pd_k,2),2);

% similarmente para outros:
fp_k = squeeze(mean(FPRt,1));
tn_k = squeeze(mean(TNRt,1));
fn_k = squeeze(mean(FNRt,1));

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
