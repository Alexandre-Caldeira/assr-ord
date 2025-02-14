function [t_decisao,TP,FP,TN,FN] = run_betamsc_cgst(...
    intens,...
    volunt,...
    M, ...
    Mstep, ...
    K, ...
    Kmax, ...
    signal_freq_bins, ...
    noise_freq_bins,...
    SIGNALS)
%run_betamsc_cgst Summary of this function goes here
%   Detailed explanation goes here
%  [t_decisao,TP,FP,TN,FN] = run_betamsc_cgst(intens, volunt, M, %disj, K,freq_bins)
% 
%       temos: M, K, alpha, gamma, intensidade, voluntario
%       calculamos: VCfutilidade, VCeficacia
%
%       aplicamos: 
%   para cada k=1:K
%     em MSC de recorte com janela de tamanho M = M + k*Mstep:
%      tki = 1:Mstep:Mlimite-Mmin; ?
%      tkf = Mmin:Mstep:Mlimite;   ?
%     testamos se: 
%         tem valor somado acumulado >= VCeficacia em signal_freq_bins   (TP)
%             entao ntp = ntp+1
%             e t_decisao = ?
%         tem valor somado acumulado <= VCfutilidade em signal_freq_bins (FN)
%             entao nfn = nfn+1
%         tem valor somado acumulado >= VCeficacia em noise_freq_bins    (FP)
%             entao nfp = nfp+1
%         tem valor somado acumulado <= VCfutilidade em noise_freq_bins  (TN)
%             entao ntn = ntn+1
%   


%% Inicializando variaveis
% isso aqui ta ruim
freq_bins = [signal_freq_bins, noise_freq_bins];
flag_noise = numel(signal_freq_bins)+1;
t_decisao = zeros(Kmax,numel(freq_bins));
TP = zeros(Kmax,numel(freq_bins));
FP = zeros(Kmax,numel(freq_bins));
TN = zeros(Kmax,numel(freq_bins));
FN = zeros(Kmax,numel(freq_bins));

%% Calcula ou busca valores criticos

% Refatorar para programacao dinamica
[aThresholds, gThresholds] = vc_beta_cgst(M, K);

%% Aplica testes

exam = zeros(numel(freq_bins),K);

if size(SIGNALS,2)>=M+K*Mstep
for k=1:K

    % Calcula indices dos recortes
    ind_final = M+k*Mstep;
    ind_inicial = ind_final-M+1;
    recorte = SIGNALS(:,ind_inicial:ind_final);
    Ps = msc_fft(recorte,M);

    for idx_freq = 1:numel(freq_bins)
    
        exam(idx_freq,k) = Ps(freq_bins(idx_freq));
    
        % SIGNAL
        if idx_freq < flag_noise && ...                         % not noise
            sum(exam(idx_freq,1:k)) >= aThresholds(k)    % detected

            TP(k,idx_freq) = TP(k,idx_freq)+1;
            t_decisao(k,idx_freq) =~sum(t_decisao(:,idx_freq));
    
        elseif idx_freq < flag_noise && ...                     % not noise
                sum(exam(idx_freq,1:k)) <= gThresholds(k)% gave up 

            FN(k,idx_freq) = FN(k,idx_freq)+1;
            t_decisao(k,idx_freq) =~sum(t_decisao(:,idx_freq));

        % NOISE
        elseif idx_freq >= flag_noise && ...                    % is noise
                sum(exam(idx_freq,1:k)) >= aThresholds(k)       % detected

            FP(k,idx_freq) = FP(k,idx_freq)+1;
            t_decisao(k,idx_freq) =~sum(t_decisao(:,idx_freq));
            
        elseif idx_freq >= flag_noise && ...                    % is noise
                sum(exam(idx_freq,1:k)) <= gThresholds(k)       % gave up

            TN(k,idx_freq) = TN(k,idx_freq)+1;
            t_decisao(k,idx_freq) = ~sum(t_decisao(:,idx_freq));
            
        end

    end

end
end
%% Diagnostico ou retorno



end