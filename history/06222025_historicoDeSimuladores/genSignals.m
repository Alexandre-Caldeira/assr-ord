function [S1, S2, S3, S4, S5, S6] = genSignals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas)
% GENSIGNALS Função que retorna três populações de sinais no dominio da frequencia
% Saidas:
%   S1 - Sinal de ruído puro, representa o nível basal para o TFG
%   S2 - Sinal de ruído puro não janelado, representa o nível basal para o TFG
%   S3 - Sinal de ruído puro janelado, utilizado para calculo de FP
%   S4 - Sinal de ruído puro não janelado, utilizado para calculo de FP
%   S5 - Sinal janelado com determinada SNR, utilizado para calculo da taxa de detecção
%   S6 - Sinal não janelado com determinada SNR, utilizado para calculo da taxa de detecção
%
% Entradas:
%     SNRfun - Função que retorna a SNR para um sinal. (Pode ser uma constante @()1, ou uma função que retorne uma amostra de uma distribuição de probabilidades)
%         FS - Frequência de amostragem. Recomenda-se que FS seja multiplo inteiro de SFREQ 
%      SFREQ - Frequência do sinal adicionado ao ruído.
%       NFFT - Número de pontos da FFT. (Usualmente igual à FS)
%    Nsinais - Numero total de sinais 
%   Njanelas - Numero de janelas de 1 segundo por sinal
%
% Autor: Quenaz Bezerra Soares
%  Data: 25/09/2019

    NpontosTotal = NFFT*Njanelas;    % Numero total de pontos de cada sinal
    tempo = (0:NpontosTotal-1)/FS;      % Vetor de tempo utilizado para gerar o sinal    
    
    % Pre-aloca as matrizes
    S1 = zeros(NFFT/2,Njanelas,Nsinais);
    S2 = zeros(NpontosTotal/2,1,Nsinais);
    S3 = zeros(NFFT/2,Njanelas,Nsinais);
    S4 = zeros(NpontosTotal/2,1,Nsinais);
    S5 = zeros(NFFT/2,Njanelas,Nsinais);
    S6 = zeros(NpontosTotal/2,1,Nsinais);
    
    for ii = 1:Nsinais
        snr = SNRfun(); % Obtém a SNR para o sinal
        sigma_n = 2/NFFT;
        snr = 10^(snr/10);
        
        % Constante a ser multiplicada ao ruido e ao sinal para configurar a relação sinal ruido desejada
        SNRs = sqrt(4*sigma_n*snr/NFFT);
        SNRn = sqrt(sigma_n);
        
        % Cria o ruído que é utilizado como base para os sinais S2, S3, S4 e S5
        ruido = randn(1,NpontosTotal);  % Gera um ruido gaussiano, teoricamente de variancia unitaria e media nula
        ruido = ruido-mean(ruido);      % Força a media nula
        ruido = ruido/std(ruido)*SNRn;  % Força a variancia desejada para o sinal
        sinal = SNRs*sin(2*pi*SFREQ*tempo+rand()*2*pi)+ruido;    % Cria o sinal e adiciona o ruido
        % Normaliza o sinal e o ruído para variância unitária
        sinal = sinal./std(sinal);
        ruido = ruido./std(ruido);
        % Obtem a FFT dos sinais não janelados
        Y = fft(ruido);
        S4(:,1,ii) = Y(1,2:floor(end/2)+1);
        Y = fft(sinal);
        S6(:,1,ii) = Y(1,2:floor(end/2)+1);
        % Reorganiza a matriz para os sinais janelados
        sinal = reshape(sinal, NFFT, Njanelas);          
        ruido = reshape(ruido, NFFT, Njanelas);
        Y = fft(ruido);
        S3(:,:,ii) = Y(2:floor(end/2)+1,:);
        Y = fft(sinal);
        S5(:,:,ii) = Y(2:floor(end/2)+1,:);

        % Cria um ruído diferente S1, pois é necessário considerar que o basal foi capturado antes do sinal de teste
        ruido = randn(1,NpontosTotal);  % Gera um ruido gaussiano, teoricamente de variancia unitaria e media nula
        ruido = ruido-mean(ruido);      % Força a media nula
        ruido = ruido/std(ruido)*SNRn;  % Força a variancia desejada para o sinal
        ruido = ruido./std(ruido);
        Y = fft(ruido);
        S2(:,1,ii) = Y(1,2:floor(end/2)+1);
        ruido = reshape(ruido, NFFT, Njanelas);
        Y = fft(ruido);
        S1(:,:,ii) = Y(2:floor(end/2)+1,:);
    end
    
end