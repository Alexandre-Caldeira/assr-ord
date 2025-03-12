classdef ORDCalculator
    properties
        dataloader

        MSC
        groupMSC
        latestMSC
        multiParamMSC

        % Calculator Parameters
        startWindow
        windowStepSize
        lastWindow
        
        epochs = []
        epochs_index_metadata
        epochs_method
        K_stages    % total number of tests to be applied on the ORD
        nWindows    % number of Windows (may match K, or not)

        subjectIndices = [1:11];
        stimulusIndices = [1:5];

        startWindows = [1];
        windowStepSizes = [24 32];
        lastWindows = [50];

        % Utils
        timer   % obj timer
        id      % obj id
    end

    methods
        function obj = ORDCalculator(dataloader)
            obj.timer = tic;
            obj.id = [num2str(keyHash(obj.timer))];
            obj.dataloader = dataloader;
            
            obj.startWindow = 1;
            obj.windowStepSize = 24;
            obj.lastWindow = dataloader.duration;
            % obj.epochs = obj.startWindow:obj.windowStepSize:obj.lastWindow;

        end

        function obj = fit_epochs(obj, p)
             arguments
                obj % The ORDCalculator class

                % p: additional parameters, passed as Name-Value arguments
                %    declared below, including their default values.
                p.dataloader = obj.dataloader;
                p.channels = obj.dataloader.channels;
                p.nChannels = numel(obj.dataloader.channels);

                p.subjectIndices = obj.subjectIndices;
                p.stimulusIndices = obj.stimulusIndices;

                p.startWindows = obj.startWindows;
                p.windowStepSizes = obj.windowStepSizes;
                p.lastWindows = obj.lastWindows;

                % How will the index of the last window be computed?
                p.lastWindowCalcMethod {mustBeMember(p.lastWindowCalcMethod,...
                                    {'maxFromStart', 'maxFromLast', 'exactK', 'fromSizeType'})} = 'maxFromStart'; 
                
                % How will the number of samples in each window be computed?
                p.sizeType {mustBeMember(p.sizeType,...
                                    {'fixedSize', 'minToMax', 'minToFix', 'withResampling'})} = 'fixedSize';
                
                % Are epochs for a single ORD or bulk ORDs?
                p.single_or_bulk {mustBeMember(p.single_or_bulk,{'single', 'bulk'})} = 'single'
                p.K_stages
                p.nWindows
             end

             % Make single to bulk switch automatic (may reconsider later,
             % bad practice...!) TODO
             if numel(p.startWindows) > 1
                 p.single_or_bulk = 'bulk';
             end

             obj.subjectIndices = p.subjectIndices;
             obj.stimulusIndices = p.stimulusIndices;

             method = [p.single_or_bulk,'_',p.sizeType,'_',p.lastWindowCalcMethod];

             switch method
                 case 'bulk_fixedSize_maxFromStart' % zanoteli bulk
                     % User expects all windows to be equally sized. 
                     % User input is subjects, stimuli, stepSizes and startWindows
                     % other params will be ignored

                     % This is not quite right cell indexing allocation but
                     % works on MATLAB... fix later! TODO
                     obj.epochs = cell(numel(p.stimulusIndices), ...
                                        numel(p.subjectIndices), ...
                                        numel(p.startWindows), ...
                                        numel(p.windowStepSizes) ...
                                        );
                     obj.epochs_index_metadata = {'stimulusIndices',...
                        'subjectIndices', 'startWindows', 'windowStepSizes'};
                     obj.epochs_method = method;
                     

                     obj.nWindows = obj.epochs;
                     obj.K_stages = obj.epochs;

                     for stimulus_index = p.stimulusIndices
                         for subject_index = p.subjectIndices
                             this_exam = cell2mat(p.dataloader.groupSignals(stimulus_index,subject_index));
                             this_lastWindow = size(this_exam, 2);

                             for this_startWindow = p.startWindows
                                 for this_windowStepSize = p.windowStepSizes

                                     this_epoch = this_startWindow:this_windowStepSize:this_lastWindow;
                                     
                                     obj.epochs{stimulus_index, ...
                                         subject_index, ...
                                         this_startWindow, ...
                                         this_windowStepSize} = this_epoch;
                                     
                                     obj.nWindows{stimulus_index, ...
                                         subject_index, ...
                                         this_startWindow, ...
                                         this_windowStepSize} = numel(this_epoch);

                                     obj.K_stages{stimulus_index, ...
                                         subject_index, ...
                                         this_startWindow, ...
                                         this_windowStepSize} = obj.nWindows(stimulus_index, ...
                                                                             subject_index, ...
                                                                             this_startWindow, ...
                                                                             this_windowStepSize);
                                end
                             end

                         end
                     end

                 case 'bulk_fixedSize_maxFromLast' % inverted zanoteli
                     % User expects all windows to be equally sized. 
                     % User input is subjects, stimuli, stepSizes and lastWindows
                     % User may pass startWindows
                     % other params will be ignored

                     obj.epochs = cell(numel(p.stimulusIndices), ...
                                        numel(p.subjectIndices), ...
                                        numel(p.startWindows), ...
                                        numel(p.windowStepSizes), ...
                                        numel(p.lastWindows)...
                                        );
                     obj.epochs_index_metadata = {'subjectIndices',...
                        'stimulusIndices', 'startWindows', ...
                        'windowStepSizes', 'lastWindows'};
                     obj.epochs_method = method;

                     obj.nWindows = obj.epochs;
                     obj.K_stages = obj.epochs;

                     for stimulus_index = p.stimulusIndices
                         for subject_index = p.subjectIndices
                             for this_firstWindow = p.startWindows
                                 for this_windowStepSize = p.windowStepSizes
                                    for this_lastWindow = p.lastWindows
                                        this_epoch = flip(p.this_lastWindow:-p.windowStepSize:this_firstWindow);
                                 
                                         obj.epochs{stimulus_index, ...
                                             subject_index, ...
                                             this_firstWindow, ...
                                             this_windowStepSize, ...
                                             this_lastWindow} = this_epoch;
                                         
                                         obj.nWindows{stimulus_index, ...
                                             subject_index, ...
                                             this_firstWindow, ...
                                             this_windowStepSize, ...
                                             this_lastWindow} = numel(this_epoch);

                                         obj.K_stages{stimulus_index, ...
                                             subject_index, ...
                                             this_firstWindow, ...
                                             this_windowStepSize, ...
                                             this_lastWindow} = p.nWindows;
                                    end
                                 end
                             end
                         end
                     end

                 case 'bulk_fixedSize_exactK' % chesnaye bulk published
                     % Requires 1 or more first windows and Ks
                     % User expects exactly K windows, starting at defined point 
                     % User input is subjects, stimuli and K
                     % User may pass startWindows
                     % other params will be ignored

                     obj.epochs = cell(numel(p.stimulusIndices), ...
                                         numel(p.subjectIndices), ...
                                         numel(p.startWindows), ...
                                         numel(p.K_stages));
                     obj.epochs_index_metadata = {'subjectIndices',...
                        'stimulusIndices', 'startWindows','K_stages'};
                     obj.epochs_method = method;

                     obj.nWindows = obj.epochs;

                     for stimulus_index = p.stimulusIndices
                         for subject_index = p.subjectIndices
                             this_exam = cell2mat(p.dataloader.groupSignals(stimulus_index,subject_index));
                             this_lastWindow = size(this_exam, 2);
                             
                             for this_firstWindow = p.startWindows %!
                                 for this_K = p.K_stages
                                     
                                     obj.epochs{stimulus_index, ...
                                         subject_index, ...
                                         this_firstWindow, ...
                                         this_K} = ceil(linspace(this_firstWindow,this_lastWindow,this_K));
                                     
                                     obj.nWindows{stimulus_index, ...
                                         subject_index, ...
                                         this_firstWindow, ...
                                         this_K} = this_K;

                                     obj.K_stages{stimulus_index, ...
                                         subject_index, ...
                                         this_firstWindow, ...
                                         this_K} = this_K;
                                 end
                             end
                         end
                     end

                 case 'bulk_fizedSize_fromSizeType' % chesnaye bulk max tests
                     % User expects all windows to be equally sized. 
                     % User input is subjects, stimuli and stepSizes
                     % User may pass startWindows
                     % other params will be ignored
                     obj.epochs = cell(numel(p.stimulusIndices), ...
                                         numel(p.subjectIndices), ...
                                         numel(p.startWindows), ...
                                         numel(p.windowStepSizes));

                     obj.epochs_index_metadata = {'subjectIndices',...
                        'stimulusIndices', 'startWindows',...
                        'windowStepSizes'};
                     obj.epochs_method = method;

                     obj.nWindows = obj.epochs;
                     obj.K_stages = obj.epochs;

                     for stimulus_index = p.stimulusIndices
                         for subject_index = p.subjectIndices
                             this_exam = cell2mat(p.dataloader.groupSignals(stimulus_index,subject_index));    
                             this_lastWindow = size(this_exam, 2);
    
                             for this_firstWindow = p.startWindows %!    
                                 for this_windowStepSize = p.windowStepSizes

                                     this_epoch = this_firstWindow:this_windowStepSize:this_lastWindow;

                                     obj.epochs{stimulus_index, ...
                                         subject_index, ...
                                         this_firstWindow, ...
                                         this_windowStepSize} = this_epoch;
                                    
                                     obj.nWindows{stimulus_index, ...
                                         subject_index, ...
                                         this_firstWindow, ...
                                         this_windowStepSize} = numel(this_epoch);

                                     obj.K_stages{stimulus_index, ...
                                         subject_index, ...
                                         this_firstWindow, ...
                                         this_windowStepSize} = p.nWindows;
                                 end
                             end
                         end
                     end


                 case 'single_fixedSize_maxFromStart' % zanoteli single
                    obj.epochs = p.startWindows:p.windowStepSizes:p.lastWindows; 
                    obj.epochs_index_metadata = 'single-exam';
                    obj.epochs_method = method;
                    obj.nWindows = numel(obj.epochs)-1;
                    obj.K_stages = obj.nWindows;

                 case 'single_fixedSize_fromSizeType' % chesnaye single
                    obj.epochs = ceil(linspace(p.startWindows,p.lastWindows,p.K_stages));
                    obj.epochs_index_metadata = 'single-exam';
                    obj.epochs_method = method;
                    obj.nWindows = p.K_stages;
                    obj.K_stages = p.K_stages;

                 otherwise
                     error('This sizeType-lastWindowCalcMethod method pair was not implemented yet.')
             end
             
        end

        function obj = compute_msc(obj, p)
            arguments
                obj % The ORDCalculator class

                % p: additional parameters, passed as Name-Value arguments
                %    declared below, including their default values.
                p.dataloader = obj.dataloader;
                p.channels = obj.dataloader.channels;
                p.nChannels = numel(obj.dataloader.channels);

                p.startWindow = obj.startWindow;
                p.windowStepSize = obj.windowStepSize;
                p.lastWindow = obj.lastWindow;
                
                p.epochs = obj.epochs;

                % epochCalcMethod define how
                p.epochCalcMethod  {mustBeMember(p.epochCalcMethod,...
                                    {'zanoteli','chesnaye'})}  = 'zanoteli'
                p.nWindows = obj.nWindows;
                p.K_stages = obj.K_stages;
                                
            end

            if isempty(p.epochs)
                switch p.epochCalcMethod
                    case 'zanoteli'
                        % obj.epochs = p.startWindow:p.windowStepSize:p.lastWindow;
                        % p.nWindows = numel(obj.epochs)-1;
                        % p.K_stages = p.nWindows;
                        p.lastWindowCalcMethod = 'maxFromStart';
                        p.sizeType = 'fixedSize';
    
                    case 'chesnaye'    
                        % obj.epochs = ceil(linspace(p.startWindow,p.lastWindow,p.K_stages));
                        % p.nWindows = p.K_stages;
                        p.lastWindowCalcMethod = 'fromSizeType';
                        p.sizeType = 'fixedSize';
                end

                obj = obj.fit_epochs( ...
                    lastWindowCalcMethod=p.lastWindowCalcMethod, ...
                    sizeType = p.sizeType, ...
                    startWindows = obj.startWindow, ...
                    windowStepSizes = obj.windowStepSize,...
                    lastWindows = obj.lastWindow,...
                    K_stages = p.K_stages...
                    );
                
                p.nWindows = obj.nWindows;
                p.K_stages = obj.K_stages;

            end
            
            % IS THIS THE CORRECT/BEST PRACTICE ACCESS??
            Y  = p.dataloader.SIGNALS;           
            obj.MSC = zeros([p.dataloader.nBins, p.nWindows, p.nChannels]);

            for channel = p.channels
                for epoch_index = 1:p.nWindows-1
                    epochStart = p.epochs(epoch_index);
                    % if epoch_index+1 > p.nWindows
                    %     disp('aqui')
                    % end

                    epochEnd = p.epochs(epoch_index+1)-1;
                    this_windowStepSize = epochEnd - epochStart+1;

                    if epochEnd > size(Y,2)
                        disp('aqui')
                    end
                    current_epoch = squeeze(Y(:,epochStart:epochEnd,channel));
                    
                    obj = obj.zanotelli_msc_fft(current_epoch, this_windowStepSize);
                    obj.MSC(:, epoch_index, channel) = obj.latestMSC;
                end
            end

        end


        % STILL UNDER DEVELOPMENT
        function obj = bulk_compute_msc(obj,p)
            arguments
                obj % The ORDCalculator class

                % p: additional parameters, passed as Name-Value arguments
                %    declared below, including their default values.
                p.dataloader = obj.dataloader;
                p.channels = obj.dataloader.channels;
                p.nChannels = numel(obj.dataloader.channels);

                % parameters for exam reload (should this be here?)
                p.subjectIndices = obj.subjectIndices;
                p.stimulusIndices = obj.stimulusIndices; 
                p.epochs = obj.epochs;
                                
            end

            % obj.parameterizedMSC
            obj.groupMSC = cell(size(p.epochs));

           for stimulusIndex = p.stimulusIndices
                for subjectIndex = p.subjectIndices
                    selected_epochs = p.epochs(stimulusIndex, subjectIndex,:);

                    for params_idx = 1:numel(selected_epochs)
                        if ~isempty(cell2mat(selected_epochs(params_idx)))

                            current_epoch = cell2mat(selected_epochs(params_idx));

                            epoch_idx = sub2ind(size(obj.epochs), ...
                                stimulusIndex, subjectIndex, params_idx );

                            p.nWindows = cell2mat(obj.nWindows(epoch_idx));

                            p.dataloader.SIGNALS = cell2mat(p.dataloader.groupSIGNALS(stimulusIndex,subjectIndex));
                            p.dataloader.nBins = size(p.dataloader.SIGNALS, 1);

                            obj = obj.compute_msc( ...
                                    dataloader = p.dataloader,...    
                                    channels = p.channels,...
                                    nChannels = p.nChannels, ...
                                    nWindows = p.nWindows, ...
                                    epochs = current_epoch ...
                                    );

                            obj.groupMSC{epoch_idx} = obj.MSC;


                        end 
                    end
                end
           end

        end

        function obj = compute_msc_on_all_channels(obj,varargin)
            if nargin>1
                obj.windowStepSize = varargin{1};
            end

            Y  = obj.dataloader.SIGNALS;            

            % Full MSC matrix is nBins x M x nChannels
            obj.nWindows = floor((size(Y,2)-1)/obj.windowStepSize);
            obj.epochs = obj.startWindow:obj.windowStepSize:obj.lastWindow;
            obj.MSC = zeros([obj.dataloader.nBins, obj.nWindows, ...
                                    obj.dataloader.nChannels]);
            
            for channel = obj.dataloader.channels
                for window_index = 1:numel(obj.epochs)-1
                    epochStart = obj.epochs(window_index);
                    epochEnd = obj.epochs(window_index+1)-1;                     
                
                    current_epoch = squeeze(Y(:,epochStart:epochEnd,channel));
    
                    obj = obj.zanotelli_msc_fft(current_epoch, obj.windowStepSize);
                    obj.MSC(:, window_index, channel) = obj.latestMSC;
                end
            end
        end
                

        function obj = zanotelli_msc_fft(obj,Y,M)
            if (size(Y,2) ~= M)
                error('Tamanho da janela diferente'); 
            end

            % y = y(1:tamanho_janela*M,1); 
            % %dividir em janela; 
            % y =reshape(y,tamanho_janela,M); 
            % 
            % %aplicar a fft; 
            % Y =fft(reshape(y,tamanho_janela,M)); %
            
            %MSC
            obj.latestMSC =  abs(sum(Y,2)).^2./(M*sum(abs(Y).^2,2));
        end


        function age(obj)
            fprintf( ...
                '\n\tThis ORDCalculator was built %0.2f seconds ago.\n\n', ...
                round(toc(obj.timer),2))
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Not working! Datastructure bugs:
        function [KA,F,KAcrit] = felix_aMSC(obj,tj,alpha)
            y = obj.dataloader.signals;
            fs = obj.dataloader.fs;

            %
            %Averaged magnitude-squared coherence (aMSC)
            %
            %Extension of the Magnitude-squared coherence (MSC) concepts to multi
            %channel analysis by averaging together many different MSCs. This is a
            %MORD technique that makes use of both magnitude and phase and can be used
            %as a detector of hidden periodicities in noise, since N real-valued
            %signals are available. The method is described in "not yet." Note: The
            %critical values were obtained via Monte Carlo simulations (1e6 iterations)
            %for a maximum of 800 epochs and alpha=0.01, 0.05 and 0.1.
            %
            %Sintaxes:
            %
            %KA = aMSC(y,tj) => aMSC spectrum.
            %[KA,F] = aMSC(y,tj,fs) => aMSC spectrum plus frequency vector.
            %[KA,F,KAcrit] = aMSC(y,tj,fs,alpha) =>  returns also the critical value.
            %
            %Input Parameters:
            %
            %y => matrix whose columns are the signals.
            %tj => number of points of each epoch in which the signals will be divided.
            %fs => sample rate of signals.
            %alpha => significance level, e.g. alpha = 0.05.
            %
            %Example: aMSC using two signals
            %
            %y1 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-15,'measured','db')';
            %y2 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-16,'measured','db')';
            %[KA,F,KAcrit] = aMSC([y1 y2],100,100,0.05);
            %figure;plot(F,KA,'b',21,1.05*KA(22),'-rv', [0 50],[KAcrit KAcrit],'k--')
            %axis([0 F(end) 0 1]);xlabel('Frequency (Hz)');ylabel('aMSC')
            
            % From Leonardo Bonato Felix - Feb/2019
            % adapted by Alexandre Gomes Caldeira - Mar/2025
            
            
            [tamsinal,N] = size(y);
            nfft = fix(tj/2)+1;  %number of points for FFT
            M = fix(tamsinal/tj); %number of windows
            y = y(1:M*tj,:); %Prevents shorter windows
            
            %Start of aMSC algorithm ---------------
            yfft = fft(reshape(y,[],M,N));
            
            %MSC
            MSC1 = zeros(nfft,N);
            for jj = 1:N
                
                MSC1(1:nfft,jj) = abs(sum(yfft(1:nfft,1:M,jj),2)).^2./(M*sum(abs(yfft(1:nfft,1:M,jj)).^2,2));
                
            end
            
            KA = mean(MSC1,2); %aMSC
            
            if nargin > 2
                F = 0:fs/tj:fs/2; %frequency vector
                if nargin == 4
                    nRuns = 10000;
                    y  = betarnd(1,M-1,nRuns,N);
                    aux = mean(y,2);
                    KAcrit = quantile(aux,1-alpha);
                    %KAcrit = aMSCcrit(N,M,alpha); %simulated 1e6 critical value
                end
            end
        end

        % Not working! Datastructure bugs:
        function [k2N,F,k2Ncrit] = felix_MMSC(obj,tj,alpha)
            y = obj.dataloader.signals;
            fs = obj.dataloader.fs;

            %
            %Multiple magnitude-squared coherence (MMSC)  
            %
            %MORD technique that makes use of both magnitude and phase and can be used
            %as a detector of hidden periodicities in noise, since N real-valued 
            %signals are available. If only one signal is used, than it reduces to the
            %standart MSC. The method was first described in "Miranda de Sa, AMFL, Felix, LB 
            %and Infantosi, AFC (2004). A Matrix-based Algorithm for Estimating 
            %Multiple Coherence of a Periodic Signal and its Application to the 
            %Multi-channel EEG during Sensory Stimulation. IEEE Transactions on 
            %Biomedical Engineering. 51(7):1140-6." The implementation seen in this 
            %code was done by refering to "Netto AD., Infantosi AFC, Miranda de Sá 
            %AMFL(2015). A Sweep Operator-Based Algorithm for Multiple Coherence 
            %Estimation in BCI Applications. In: Lackovi? I., Vasic D. (eds) 6th 
            %European Conference of the International Federation for Medical and 
            %Biological Engineering. IFMBE Proceedings, vol 45. Springer, Cham 
            %
            %Sintaxes:
            %
            %k2N = MMSC(y,tj) => MMSC spectrum.
            %[k2N,F] = MMSC(y,tj,fs) => MMSC spectrum plus frequency vector. 
            %[k2N,F,k2Ncrit] = MMSC(y,tj,fs,alpha) =>  returns also the theoretical 
            %critical value.
            %
            %Input Parameters:
            %
            %y => matrix whose columns are the signals.  
            %tj => number of points of each epoch in which the signals will be divided.
            %fs => sample rate of signals.
            %alpha => significance level, e.g. alpha = 0.05.
            %
            %Example: MMSC using two signals 
            %
            %y1 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-15,'measured','db')'; 
            %y2 = awgn(sin(2*pi*21*(linspace(0,10,1000))),-16,'measured','db')';
            %[K2N,F,K2Ncrit] = MMSC([y1 y2],100,100,0.05);
            %figure;plot(F,K2N,'b',21,1.05*K2N(22),'-rv', [0 50],[K2Ncrit K2Ncrit],'k--')
            %axis([0 F(end) 0 1]);xlabel('Frequency (Hz)');ylabel('MMSC')
            
            % From Bonato mar/2019, 
            % adapted by Caldeira mar/2025
            
            nf = fix(tj/2)+1;
            N = size(y,2); %numero de sinais
            M = fix(size(y,1)/(tj));   %determina numero de segmentos
            Sfft = zeros(tj,N,M);   %acumulador para fft em anel
            k2N = zeros(nf,1);  %calculo da coerencia
            idl = 1;
            
            %monta matriz de espectros
            for k=1:M
                
                %retira tendência linear
                %Sfft(:,:,k) = fft(detrend(canais(idl:(idl+L-1),:),'linear'));
                Sfft(:,:,k) = fft(y(idl:(idl+tj-1),:));
                
                idl = idl+tj;
            end
            
            for kf =1:nf
                %reinicia matrizes para calculo
                %V = zeros(N,1);
                %Specm = zeros(N,N);
                
                %matriz de espectros aumentados, parte com V e VH
                Specm_a = zeros(N+1,N+1);  
                for p = 1:N
                    for ks = 1: M
                        %V(p) = V(p) + Sfft(kf,p,ks);
                        Specm_a(N+1,p) = Specm_a(N+1,p) + Sfft(kf,p,ks);
                    end
                    %VH - hermitiano de V
                    Specm_a(p,N+1) = conj(Specm_a(N+1,p));
                end
                
                
                %monta mtriz de espectro aumentado
                for p = 1:N
                    for q = 1:p
                        for ks=1:M
                            Specm_a(p,q) = Specm_a(p,q) + conj(Sfft(kf,p,ks)).*Sfft(kf,q,ks);
                        end
                        Specm_a(q,p) = conj(Specm_a(p,q));
                    end
                end
                Specm_a(N+1,N+1) = 1;
                Specm_as = Msweep(Specm_a,N);
                k2N(kf) = (1-real(Specm_as(N+1,N+1)))/M;
            end
            
            if nargin > 2
                F = (0:(nf-1))*fs/tj; %frequency vector
                if nargin == 4
                    Fcrit = finv(1-alpha,2*N,2*M-2*N);% F CDF
                    k2Ncrit = Fcrit/(((M-N)/N)+Fcrit); %critical value
                end
            end               
            
        end

        
    end
    
end