classdef PreProcessor
    % TODO: add all filters from (Antunes et al., 2016?) for further use 
    properties
        filterFcLower = 70;                  % Lower cutoff frequency (Hz)
        filterFcUpper = 1000/2 - 1;          % Upper cutoff frequency (Hz)
        filterOrder   = 8;                   % Butterworth filter order
        fs = 1000;
        

        zanoteliGain  = 200;
        zanoteliCut
        
        processedSignals
    end

    methods
        function obj = PreProcessor(dataloader)
            obj.fs = dataloader.fs;
            obj.filterFcUpper = obj.fs/2 - 1;
            obj.zanoteliCut = 0.1/obj.zanoteliGain; % for artifact removal
        end

        function obj = zanoteliPreProcessing(obj, dataloader)
            % preprocessSignal removes the DC offset and applies a Butterworth bandpass filter.
            % expected datastructure is signal = zeros(fs, epochs, nchannels)
            

            % TODO: remover ALTISSIMO acoplamento aqui! deveria ser
            % interface na entrada e nao consumir o obj inteiro ou sla
            for channel = dataloader.channels
                % If mock needed: dataloader.signal = nan(dataloader.fs, dataloader.duration, dataloader.nChannel)
                x = squeeze(dataloader.signals(:,:,channel));
    
                %1segundo de sinal 
                nfft = dataloader.nfft;
                 
                %retirar componente DC por janela 
                % (fiz isso pq no processamento em tempo real é por janela)
                % tirar a média de cada de cada trecho - devido a remoção
                x = x - repmat(mean(x),nfft,1); 
                    
                % Excluir os dois primeiros segundos do inicio da coleta 
                x(:,1:2,:) = []; 
                    
                % Encontrar o valor máximo por canal 
                Vmax = max(abs(x),[],1);

                % Remover o ruído amplitude 
                ind = Vmax > obj.zanoteliCut;
                x = x(:,~ind);  

                % Limitar o tamanho para o valor máximo. 
                Mmax = dataloader.zanoteliSuggestedMMax(dataloader.zanoteliStimulusIndex);
                obj.processedSignals(:,:,channel) = x(:,1:Mmax);

            end

            % Remove DC offset
            % obj.processedSignals = signal - mean(signal, 2);
            % 
            % [b, a] = butter(obj.filterOrder, ... 
            %     [obj.filterFcLower, obj.filterFcUpper] / (obj.fs/2)); 
            % obj.processedSignals = filtfilt(b, a, obj.processedSignals); 

        end
    end
end