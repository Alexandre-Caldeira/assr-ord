function [RES, RESESP] = prepareVP(RAIZ, M)
% Fazer testes mais simples:
%   Utilizar apenas o eletrodo FCz,30 janelas
%   50, 60 e 70 dB
% clear
% clc
% RAIZ = 'C:\Users\quena\Desktop\Dados para usar na Monografia\Dados_organizados_Colatina\';
%Njanelas = [30 60];
% M = 60;

dir(RAIZ)
res = nan(14,1,8*11);
RES = nan(14,1,8*11*16);
resESP = nan(14,1,11*20);
RESESP = nan(14,1,11*20*16);
filtros = {'*ESP.mat','*_70dB.mat'};

for canal = 1:16
    S3 = nan(500, M, 11);
    S4 = nan(500*M,1,11);
    for iEST = 1:2
        arquivos = dir([RAIZ filtros{iEST}]);
        arquivos = {arquivos.name};
        for iIND = 1:11
            dados = load([RAIZ filesep arquivos{iIND}]);
            FS = dados.Fs;
            xtemp = dados.x;
            %binsM = dados.binsM;% Estava pegando os bins errados
            binsM = dados.freqEstim;
            if iEST==1 % Tratamento especial para o EEG Espontâneo
                [~,n,~] = size(xtemp);
                xEsp = xtemp(:,randperm(n),:); % Embaralha o sinal
                temp = xEsp(:,:,canal);
                x = forceSize(temp(:),FS*M);
                x = x - mean(x);
                x = x ./ std(x);
                if FS ~=1000
                    [p,q] = rat(1000 / FS);
                    x = resample(x,p,q);
                    x = x - mean(x);
                    x = x ./ std(x);
                end
                Y = fft(x);%*sqrt(1000/FS);
                S4(:,1,iIND) = Y(2:500*M+1,:);
                x = reshape(x, 1000, M);
                Y = fft(x);%*sqrt(1000/FS);
                S3(:,1:M,iIND) = Y(2:501,:);
            end
            S1 = nan(M*500,1,1);
            S2 = nan(500,M,1);
            temp = xtemp(:,:,canal);
            x = forceSize(temp(:),FS*M);
            x = x - mean(x);
            x = x ./ std(x);
            if FS ~=1000
                [p,q] = rat(1000 / FS);
                x = resample(x,p,q);
                x = x - mean(x);
                x = x ./ std(x);
            end
            Y = fft(x);%*sqrt(1000/FS);
            S1(:,1,1) = Y(2:floor(end/2)+1);
            x = reshape(x, 1000, M);
            Y = fft(x);%*sqrt(1000/FS);
            S2(:,:,1) = Y(2:floor(end/2)+1,:);
            for dna = 1:14
                if iEST==1
                    for iBIN = 1:20
                        resESP(dna,1,(iIND-1)*20+iBIN) = funcoesPrimitivas(dna, S2, S3(:,:,iIND), S1, S4(:,:,iIND), M, iBIN+80);
                    end
                else
                    for iBIN = 1:8
                        res(dna,1,(iEST-2)*11*8+(iIND-1)*8+iBIN) = funcoesPrimitivas(dna, S2, S3(:,:,iIND), S1, S4(:,:,iIND), M, binsM(iBIN));
                    end
                end
            end
            disp((iIND/11+iEST-1)/6*100)
        end
    end
    RESESP(:,:,(canal-1)*11*20+1:canal*11*20) = resESP;
    RES(:,:,(canal-1)*8*11+1:canal*8*11) = res;
end

% Seleciona os canais
for canal = [15 13:-1:2]
    RES(:,:,(canal-1)*8*11+1:canal*8*11) = [];
end



    function x = forceSize(x,n)
        assert(sum(size(x)~=1)==1,'x deve ser um vetor linha ou coluna')
        assert(length(size(x))==2,'x deve ser um vetor linha ou coluna')
        [MM,NN]= size(x);
        dN = n-max(MM,NN);
        if dN>0
            x(end+1:end+dN) = nan(1,dN);
        elseif dN<0
            x(end+dN+1:end) = [];
        end
    end

end