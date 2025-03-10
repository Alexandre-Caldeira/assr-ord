classdef DataLoader
% IN: experimental .mat data path or MC sim parameters
% out: Y (freq. data for testing)

    properties
       
        signals % y (timeseries data, MC sim or exp aquisition)
        SIGNALS % frequency data
        signalFrequencies = [82    84    86    88    90    92    94    96];
        noiseFrequencies = [];
        
        % Processing Parameters
        channels = 1:16;  % Index for EEG leads (channels) in data
        fs       = 1000;  % Sampling frequency (Hz)
        nfft              % Number of FFT points per epoch (1 sec window in NIASv1)
        nBins             % Number of 'positive frequency' values on the spectrum
        nChannels         % length of channels vector
        duration = [60];  % exam duration [secs], per zanoteliStimulusIndex intensity
        totalSamples

        % For MC Simulation        
        noiseMean = -10
        noiseStd = 1
        SNRfun

        % For ASSR EEG Data
        zanoteliSubjectIndex = 1;
        zanoteliStimulusIndex = 1;

        zanoteliSuggestedMMax = [50; 50; 240; 440; 440; 20];
        zanoteliStimulusNames = {'70dB'; '60dB'; '50dB';
                                 '40dB'; '30dB';'ESP'};
        zanoteliSubjects  = {'Ab'; 'An'; 'Bb'; 'Er'; 'Lu';...
                             'So'; 'Qu'; 'Vi'; 'Sa'; 'Ti';'Wr'};
        zanoteliLeads = {'FC'; 'F4'; 'T6'; 'P4'; 'T4'; 'Oz'; 'C4'; 'T5';...
                         'P3'; 'F7'; 'F3'; 'T3'; 'C3'; 'Fz'; 'Pz'; 'Cz'}

        % Utils
        timer   % obj timer
        id      % obj id
        
        inPath  = 'C:\PPGEE\SBEB_CBA_24\CGST_figuras\Sinais_EEG\' % where y comes from
        outPath % where Y is saved
        mode    % simulation mode 'sim' or 'exp'
        

    end
    
    methods
        function obj = DataLoader(mode,varargin)
            obj.timer = tic;
            obj.id = [num2str(keyHash(keyHash(mode)+obj.timer))];
            obj.outPath = [pwd,'\'];
            
            if nargin == 2
                obj.inPath = varargin{1};
            elseif nargin == 3
                obj.inPath = varargin{1};
                obj.outPath = varargin{2};
            end
    
            obj.nfft = obj.fs; 
            obj.nBins = floor(obj.fs/2)+1; % ...
            obj.noiseFrequencies = randi([300 500],1, ...
                                    numel(obj.signalFrequencies));

            % By default, SNR is:
            % Adiciona ruido gaussiano branco com SNR aleatoria 
            obj.SNRfun = @() obj.noiseMean+obj.noiseStd()*randn(1);

            obj.mode = mode;
            if matches(obj.mode(1:3),"sim", IgnoreCase=true)
                
                % Run MC Sim
                obj = obj.genSimulatedSignals();
                obj = obj.computeFFT();

            elseif matches(obj.mode(1:3),"exp", IgnoreCase=true)
                % Try to load exp data for validation

                %carregar o voluntário
                zanoteliSubjectIndex = obj.zanoteliSubjectIndex; 
                
                %intensidadde 
                zanoteliStimulusNames = cell2mat(obj.zanoteliStimulusNames(obj.zanoteliStimulusIndex));
                
                obj = obj.loadEEGData(zanoteliSubjectIndex, zanoteliStimulusNames);

            else
                % Throw error
                error('DataLoader mode input is invalid.'); 
            end

        end

        function obj = loadEEGData(obj, zanoteliSubjectIndex, zanoteliStimulusNames)
            % Check if obj is properly constructed
            if isempty(obj.inPath) 
                error('EEG data file path not specified.'); 
            end 
            
            subject_id = cell2mat(obj.zanoteliSubjects(zanoteliSubjectIndex));

            % Build filepath and load
            % data = load(filepath); 
            data = load([obj.inPath subject_id zanoteliStimulusNames], ...
                'x','Fs','binsM','freqEstim');
            
            if isfield(data, 'x') 
                obj.signals = data.x; 
                obj.totalSamples = size(obj.signals,2)*size(obj.signals,1);
            else 
                error('The EEG data file must contain variable "x".'); 
            end 

            if isfield(data, 'Fs') 
                obj.fs = data.Fs;
                obj.nfft = obj.fs;
                obj.nBins = floor(obj.fs/2)+1;
            else 
                error('The EEG data file must contain variable "x".'); 
            end

        end

        function inspectExam(obj)
            fprintf(...
            '\n\tExam is %s stimulus on subject %s,\n\t measuring on %s for %s seconds.\n\n', ...
                cell2mat(obj.zanoteliStimulusNames(obj.zanoteliStimulusIndex)), ...
                cell2mat(obj.zanoteliSubjects(obj.zanoteliSubjectIndex,:)),...
                cell2mat(obj.zanoteliLeads(1)),...
                num2str(size(obj.signals, 2))...
                )
        end

        function obj =resetSNRfun(obj,noiseMean,noiseStd)
            % Alguns exemplos:
            % noise_var_mean = 2^2;
            % noise_var = @() (noise_var_mean+ sqrt(noise_var_var)*randn(1));
            % noise_sd = @() randi([1,45],1)/10;

            % Test case:
            % dtl = DataLoader('sim');
            % vec = zeros(1,1e5);
            % for i = 1:numel(vec)
            % vec(i) = dtl.SNRfun();
            % end
            % histogram(vec)

            obj.noiseMean = noiseMean;
            obj.noiseStd = noiseStd;
            obj.SNRfun = @() obj.noiseMean+obj.noiseStd()*randn(1);

        end

        function obj = genSimulatedSignals(obj) 
            % genSimulatedSignals generates simulated EEG signals (noise + sinusoidal signal). 
            % Example test usage:
            % dtl = DataLoader('sim');
            % stem(abs(dtl.SIGNALS(:,1,1)))
            % hold on
            % stem(dtl.signalFrequencies,abs(dtl.SIGNALS(dtl.signalFrequencies,1,1)),...
            % 'MarkerFaceColor','red',...
            % 'MarkerEdgeColor','red')

            obj.nChannels = numel(obj.channels);
            obj.totalSamples = obj.fs * obj.duration * obj.nChannels; 
            t = (0:obj.totalSamples-1) / obj.fs;     
        
            obj.signals = 0;
            for i = 1:numel(obj.signalFrequencies)
                fo = obj.signalFrequencies(i)-1; 
                % fo subtracted by one, such that FFT bin for frequency 
                % will match the frequency because MATLAB index starts at 1

                obj.signals = obj.signals+sin(2*pi*fo*t);
            end
           
            obj.signals = reshape(obj.signals, [obj.fs, obj.duration, obj.nChannels]);  
            for channel = 1:obj.nChannels
                for epoch = 1:obj.duration

                    x = awgn( ... % add white gaussian noise
                            obj.signals(:,epoch,channel), ... % to this section of signal
                            obj.SNRfun(),'measured','db')'; % with this SNR in measured dB

                    % Scale to uV
                    obj.signals(:,epoch,channel)  = (10^-6) * x;  
                end
            end

        end

        function obj = computeFFT(obj) 
            % computeFFT computes the FFT along each window (row-wise) and 
            % returns only the positive frequencies. 
            
            obj.duration = size(obj.signals,2);
            obj.nChannels= numel(obj.channels);
            obj.SIGNALS = zeros(obj.nBins, obj.duration, obj.nChannels);

            for channel = obj.channels
                for epoch = 1:obj.duration
                    temp = fft(squeeze(obj.signals(:,epoch,channel)),obj.nfft, 1); 
                    obj.SIGNALS(:,epoch,channel) = temp(1:obj.nBins); 
                end
            end

        end

        function obj = resetExam(obj, subjectIndex,stimulusIndex)
            obj.zanoteliSubjectIndex = subjectIndex;
            obj.zanoteliStimulusIndex = stimulusIndex;
            obj = obj.loadEEGData(obj.zanoteliSubjectIndex, ...
                         cell2mat(obj.zanoteliStimulusNames(obj.zanoteliStimulusIndex)));
            obj= obj.computeFFT();

        end

        function obj = resetSubject(obj, subjectIndex)
            % TODO: add inspect exam and recompute as optional
            obj.zanoteliSubjectIndex = subjectIndex;
            obj = obj.loadEEGData(obj.zanoteliSubjectIndex, ...
                         cell2mat(obj.zanoteliStimulusNames(obj.zanoteliStimulusIndex)));
            obj= obj.computeFFT();

        end

        function obj = resetStimulus(obj,stimulusIndex)
            obj.zanoteliStimulusIndex = stimulusIndex;
            obj = obj.loadEEGData(obj.zanoteliSubjectIndex, ...
                         cell2mat(obj.zanoteliStimulusNames(obj.zanoteliStimulusIndex)));
            obj = obj.computeFFT();

        end

        function obj = resetChannels(obj,channelsIndices)
            obj.channels = channelsIndices;
            obj.nChannels = numel(channelsIndices);
            obj = obj.computeFFT();

        end

        function obj = resetDuration(obj,newDuration)
            if obj.zanoteliSuggestedMMax(obj.zanoteliStimulusIndex) < newDuration 
                warning(['Desired duration is larger than Zanoteli suggestion.' ...
                    'This might impact on some subjects lacking data on later epochs.' ...
                    'Consider lowering exam duration.'])
            end

            obj.duration = newDuration;
            obj = obj.computeFFT();
            
        end

        function age(obj)
            fprintf('\n\tThis DataLoader was built %0.2f seconds ago.\n\n', round(toc(obj.timer),2))
        end

    end
        
end
