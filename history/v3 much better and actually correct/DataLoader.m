classdef DataLoader < handle
    % Loads or simulates EEG time-series data for a single exam configuration.
    % Manages data paths, simulation parameters, and metadata.

    properties (SetAccess = protected)
        % Core Data (for a single loaded/simulated exam)
        rawSignals % Stores time-series y(sample, window, channel)

        % Signal Parameters
        signalFrequencies = [82 84 86 88 90 92 94 96]; % Frequencies expected to contain signal
        noiseFrequencies = []; % Frequencies considered as noise (randomly selected if empty)
        fs = 1000;  % Sampling frequency (Hz)
        channels = 1:16;  % Default EEG channel indices
        nChannels % Number of channels (derived)

        % Simulation Parameters (if mode is 'sim')
        noiseMean = -10; % Mean SNR (dB) for AWGN
        noiseStd = 1;    % Standard deviation of SNR (dB) for AWGN
        simDuration = 60; % Default simulation duration (seconds)

        % Experimental Data Parameters (if mode is 'exp')
        zanoteliSubjectList = {'Ab'; 'An'; 'Bb'; 'Er'; 'Lu'; 'So'; 'Qu'; 'Vi'; 'Sa'; 'Ti'; 'Wr'};
        zanoteliStimulusNames = {'70dB'; '60dB'; '50dB'; '40dB'; '30dB'; 'ESP'};
        zanoteliLeads = {'FC'; 'F4'; 'T6'; 'P4'; 'T4'; 'Oz'; 'C4'; 'T5'; 'P3'; 'F7'; 'F3'; 'T3'; 'C3'; 'Fz'; 'Pz'; 'Cz'};
        zanoteliSuggestedMMax = [50; 50; 240; 440; 440; 20]; % Max suggested windows per stimulus

        % State
        mode % 'exp' or 'sim'
        currentSubjectIndex = 1
        currentStimulusIndex = 1
        currentSubjectName = ''
        currentStimulusName = ''
        currentDuration % Actual duration of the loaded/generated signal in windows/seconds
        isDataLoaded = false

        % Paths
        expDataPath = 'C:\PPGEE\SBEB_CBA_24\CGST_figuras\Sinais_EEG\' % Default path for experimental data
    end

    properties (Dependent)
        nSignalFrequencies
        nNoiseFrequencies
        allTestFrequencies % Combined and sorted signal and noise frequencies
    end

    methods
        function obj = DataLoader(mode, options)
            arguments
                mode {mustBeMember(mode, {'exp', 'sim'})}
                options.dataPath string = 'C:\PPGEE\SBEB_CBA_24\CGST_figuras\Sinais_EEG\'
                options.channels (1,:) {mustBeInteger, mustBePositive} = 1:16
                options.signalFrequencies (1,:) {mustBeNumeric, mustBePositive} = [82 84 86 88 90 92 94 96]
                options.noiseFrequencies (1,:) {mustBeNumeric, mustBePositive} = [] % Allow empty
                options.simDuration (1,1) {mustBeNumeric, mustBePositive, mustBeInteger} = 60
                options.noiseMean (1,1) {mustBeNumeric} = -10
                options.noiseStd (1,1) {mustBeNumeric, mustBeNonnegative} = 1
            end

            if numel(options.signalFrequencies) ~= numel(unique(options.signalFrequencies))
                error('DataLoader:InvalidInput', 'Input signalFrequencies must contain unique values.');
            end
            if ~isempty(options.noiseFrequencies) && numel(options.noiseFrequencies) ~= numel(unique(options.noiseFrequencies))
                 error('DataLoader:InvalidInput', 'Input noiseFrequencies must contain unique values if provided.');
            end

            obj.mode = mode;
            obj.expDataPath = options.dataPath;
            obj.channels = options.channels;
            obj.signalFrequencies = sort(options.signalFrequencies);
            obj.noiseFrequencies = sort(options.noiseFrequencies);
            obj.simDuration = options.simDuration;
            obj.noiseMean = options.noiseMean;
            obj.noiseStd = options.noiseStd;
            obj.nChannels = numel(obj.channels);

            if isempty(obj.noiseFrequencies)
                obj.generateNoiseFrequencies();
            end
            % Reduced verbosity
            % fprintf('[%s] DataLoader initialized in "%s" mode.\n', datetime, obj.mode);
        end

        function generateNoiseFrequencies(obj, numNoiseFreqs, range)
            arguments
                obj
                numNoiseFreqs = numel(obj.signalFrequencies)
                range = [300, 500]
            end
             if range(2) <= range(1) || range(1) < max(obj.signalFrequencies)
                warning('Invalid noise frequency range. Using default [300, 500] Hz.');
                range = [300, 500];
            end
            availableFreqs = range(2) - range(1);
            potentialNoise = range(1):range(2);
            potentialNoise = setdiff(potentialNoise, obj.signalFrequencies);

            if numNoiseFreqs > numel(potentialNoise)
                 warning('Could not find enough unique noise frequencies outside signal range. Using replacement within range.');
                 obj.noiseFrequencies = sort(randi(range, 1, numNoiseFreqs));
            else
                 obj.noiseFrequencies = sort(randsample(potentialNoise, numNoiseFreqs));
            end
             % Reduced verbosity
             % fprintf('[%s] Generated %d noise frequencies in range [%d, %d] Hz.\n', ...
             %       datetime, numel(obj.noiseFrequencies), range(1), range(2));
        end

        function success = loadExperiment(obj, subjectIndex, stimulusIndex)
            arguments
                obj
                subjectIndex {mustBeInteger, mustBePositive}
                stimulusIndex {mustBeInteger, mustBePositive}
            end
            success = false;
            if ~strcmp(obj.mode, 'exp')
                error('DataLoader is not in experimental mode.');
            end
            if subjectIndex > numel(obj.zanoteliSubjectList) || stimulusIndex > numel(obj.zanoteliStimulusNames)
                warning('Invalid subject or stimulus index.');
                return;
            end

            obj.currentSubjectIndex = subjectIndex;
            obj.currentStimulusIndex = stimulusIndex;
            obj.currentSubjectName = obj.zanoteliSubjectList{subjectIndex};
            obj.currentStimulusName = obj.zanoteliStimulusNames{stimulusIndex};
            filepath = fullfile(obj.expDataPath, [obj.currentSubjectName, obj.currentStimulusName, '.mat']);

            try
                % Reduced verbosity
                % fprintf('[%s] Loading: %s...\n', datetime, filepath);
                data = load(filepath, 'x', 'Fs');
                if isfield(data, 'x') && isfield(data, 'Fs')
                    obj.rawSignals = data.x;
                    obj.fs = data.Fs;
                    if max(obj.channels) > size(obj.rawSignals, 3)
                        warning('Requested channels (%s) exceed available channels (%d) in data. Using subset: %s.', ...
                                mat2str(obj.channels), size(obj.rawSignals, 3), mat2str(1:size(obj.rawSignals, 3)));
                        obj.channels = 1:size(obj.rawSignals, 3);
                        obj.nChannels = numel(obj.channels);
                    end
                    obj.rawSignals = obj.rawSignals(:, :, obj.channels);
                    obj.currentDuration = size(obj.rawSignals, 2);
                    obj.isDataLoaded = true;
                    success = true;
                    % Reduced verbosity
                    % fprintf('[%s] Loaded Subject %d (%s), Stimulus %d (%s). Duration: %d s, Fs: %d Hz, Channels: %d.\n', ...
                    %        datetime, subjectIndex, obj.currentSubjectName, stimulusIndex, obj.currentStimulusName, ...
                    %        obj.currentDuration, obj.fs, obj.nChannels);
                else
                    warning('Loaded file %s is missing required variables "x" or "Fs".', filepath);
                    obj.isDataLoaded = false;
                end
            catch ME
                warning('Failed to load file %s: %s', filepath, ME.message);
                obj.isDataLoaded = false;
            end
        end

        function success = generateSimulation(obj, options)
            arguments
                obj
                options.duration (1,1) {mustBePositive, mustBeInteger} = obj.simDuration
                options.noiseMean (1,1) {mustBeNumeric} = obj.noiseMean
                options.noiseStd (1,1) {mustBeNumeric, mustBeNonnegative} = obj.noiseStd
            end
            success = false;
            if ~strcmp(obj.mode, 'sim')
                error('DataLoader is not in simulation mode.');
            end

            obj.simDuration = options.duration;
            obj.noiseMean = options.noiseMean;
            obj.noiseStd = options.noiseStd;
            obj.currentDuration = obj.simDuration;
            obj.currentSubjectName = sprintf('Sim_SNR_%.1f_%.1f', obj.noiseMean, obj.noiseStd);
            obj.currentStimulusName = sprintf('%ds', obj.simDuration);
            obj.currentSubjectIndex = NaN;
            obj.currentStimulusIndex = NaN;

            % Reduced verbosity
            % fprintf('[%s] Generating simulation: Duration=%d s, SNR Mean=%.1f dB, SNR Std=%.1f dB...\n', ...
            %        datetime, obj.simDuration, obj.noiseMean, obj.noiseStd);

            totalSamplesPerChannel = obj.fs * obj.simDuration;
            t = (0:totalSamplesPerChannel-1)' / obj.fs;

            baseSignal = zeros(totalSamplesPerChannel, 1);
            for i = 1:numel(obj.signalFrequencies)
                fo = obj.signalFrequencies(i);
                baseSignal = baseSignal + sin(2*pi*fo*t);
            end

            obj.rawSignals = zeros(obj.fs, obj.simDuration, obj.nChannels);
            snrFun = @() obj.noiseMean + obj.noiseStd * randn(1);

            for ch = 1:obj.nChannels
                channelSignal = reshape(baseSignal, [obj.fs, obj.simDuration]);
                noisyChannelSignal = zeros(size(channelSignal));
                for epoch = 1:obj.simDuration
                    epochSignal = channelSignal(:, epoch);
                    currentSNR = snrFun();
                    noisyEpoch = awgn(epochSignal, currentSNR, 'measured', 'db');
                    noisyChannelSignal(:, epoch) = (1e-6) * noisyEpoch;
                end
                 obj.rawSignals(:, :, ch) = noisyChannelSignal;
            end

            obj.isDataLoaded = true;
            success = true;
            % Reduced verbosity
            % fprintf('[%s] Simulation generated.\n', datetime);
        end

        function signals = getRawSignals(obj)
            if ~obj.isDataLoaded
                error('No data loaded or generated. Call loadExperiment or generateSimulation first.');
            end
            signals = obj.rawSignals;
        end

        function clearData(obj)
            obj.rawSignals = [];
            obj.isDataLoaded = false;
            % Reduced verbosity
            % fprintf('[%s] Raw signal data cleared from DataLoader.\n', datetime);
        end

        % --- Dependent Properties ---
        function val = get.nSignalFrequencies(obj)
            val = numel(obj.signalFrequencies);
        end
        function val = get.nNoiseFrequencies(obj)
            val = numel(obj.noiseFrequencies);
        end
        function val = get.allTestFrequencies(obj)
            val = sort([obj.signalFrequencies, obj.noiseFrequencies]);
        end

        % --- Utility Methods ---
         function inspect(obj)
            fprintf('\n--- DataLoader Inspection (%s) ---\n', datetime);
            fprintf('Mode: %s\n', obj.mode);
            fprintf('Sampling Frequency (Fs): %d Hz\n', obj.fs);
            fprintf('Channels (%d): %s\n', obj.nChannels, mat2str(obj.channels));
            fprintf('Signal Frequencies (%d): %s Hz\n', obj.nSignalFrequencies, mat2str(obj.signalFrequencies));
            fprintf('Noise Frequencies (%d): %s Hz\n', obj.nNoiseFrequencies, mat2str(obj.noiseFrequencies));
            if strcmp(obj.mode, 'exp')
                fprintf('Exp. Data Path: %s\n', obj.expDataPath);
                 if obj.isDataLoaded, fprintf('Current Subject: #%d (%s), Stimulus: #%d (%s), Duration: %d s\n', obj.currentSubjectIndex, obj.currentSubjectName, obj.currentStimulusIndex, obj.currentStimulusName, obj.currentDuration); else, fprintf('Current Subject/Stimulus: None Loaded\n'); end
            else
                fprintf('Sim. Default Duration: %d s, Noise Mean: %.1f dB, Noise Std: %.1f dB\n', obj.simDuration, obj.noiseMean, obj.noiseStd);
                 if obj.isDataLoaded, fprintf('Current Simulation: Duration=%d s, SNR Mean=%.1f, SNR Std=%.1f\n', obj.currentDuration, obj.noiseMean, obj.noiseStd); else, fprintf('Current Simulation: None Generated\n'); end
            end
            fprintf('Data Loaded/Generated: %s\n', mat2str(obj.isDataLoaded));
            fprintf('------------------------------------\n');
        end
    end
end