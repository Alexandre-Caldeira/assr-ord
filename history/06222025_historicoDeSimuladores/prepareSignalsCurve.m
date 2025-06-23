NsinaisTotal = 500; %Número de sinais totais que serão gerados. Obrigatóriamente múltiplo de 100
Njanelas = [60 30];%[30 60 90 120];
NFFT = 1000;
SFREQ = 80;
FS = 1000;
Nfun = 14;
%SNR = -50:10;
SNR = -40:2:0;
nSNR = length(SNR);
SNRfun = @()-15+5*randn;    % SNR aleatória, centrada em -15, com desvio padrão igual a 5
tic
% Prepara os sinais para os números de janelas dejesados
for iJAN = 1:length(Njanelas)
    M = Njanelas(iJAN);
    % Prealoca a matriz final para o número de janelas atual
    res = nan(Nfun,nSNR,NsinaisTotal/100,100);
    resLim = nan(Nfun,nSNR,NsinaisTotal/100,100);
	for iSNR = 1:nSNR
        toc
        tic
		snr = SNR(iSNR);
		SNRfun = @()snr;
        fprintf('SNR %08.3f: ', snr)
		% Gera os sinais em grupos de 100, para evitar problemas de memória
		for iSinais = 1:NsinaisTotal/100
			[S1, S2, S3, S4, S5, S6] = genSignals(SNRfun, FS, SFREQ, NFFT, 100, M);
			for func = 1:Nfun
				res(func, iSNR, iSinais,:) = funcoesPrimitivas(func, S5, S1, S6, S2, M, SFREQ);
				resLim(func, iSNR, iSinais,:) = funcoesPrimitivas(func, S3, S1, S4, S2, M, SFREQ);
            end
            if mod(iSinais,10)==0
                fprintf('%4d janelas: %6d de %6d concluídos.\n', M, iSinais, NsinaisTotal/100);
            end
		end
	end
    RES = nan(Nfun,nSNR,NsinaisTotal);
    RESLim = nan(Nfun,nSNR,NsinaisTotal);
    for i = 1:NsinaisTotal/100
        RES(:,:,(i-1)*100+1:i*100) = res(:,:,i,:);
        RESLim(:,:,(i-1)*100+1:i*100) = resLim(:,:,i,:);
    end
    res = RES;
    save(sprintf('..\\MC_curva\\resY_%d.mat', M), 'res')
    resLim = RESLim;
    save(sprintf('..\\MC_curva\\resX_%d.mat', M), 'resLim')
end