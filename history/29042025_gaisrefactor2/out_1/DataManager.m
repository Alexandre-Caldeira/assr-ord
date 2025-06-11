% --- FILE START: DataManager.m ---
classdef DataManager
    % Manages loading of experimental data or generation of simulated data.
    % Stores metadata related to the data source.

    properties
        Mode            % 'sim' or 'exp': Specifies the data source type.
        DataSourcePath  % Path for experimental data files or simulation config.
        SavePath        % Path to save processed data or results.

        % --- Experiment Specific Properties ---
        ExpSubjects = {'Ab'; 'An'; 'Bb'; 'Er'; 'Lu'; 'So'; 'Qu'; 'Vi'; 'Sa'; 'Ti'; 'Wr'}; % Default subject list.
        ExpStimuli = {'70dB'; '60dB'; '50dB'; '40dB'; '30dB'; 'ESP'}; % Default stimulus list.
        ExpLeads = {'FC'; 'F4'; 'T6'; 'P4'; 'T4'; 'Oz'; 'C4'; 'T5'; 'P3'; 'F7'; 'F3'; 'T3'; 'C3'; 'Fz'; 'Pz'; 'Cz'}; % Default lead names.
        ExpSuggestedMMax = [50; 50; 240; 440; 440; 20]; % Max windows suggested by Zanoteli per stimulus.

        % --- Simulation Specific Properties ---
        SimSignalFrequencies = [82 84 86 88 90 92 94 96]; % Freqs (Hz) for simulated signal components.
        SimNoiseFrequencies = []; % Freqs (Hz) for noise components (can be auto-generated).
        SimNoiseMean = -10; % Default SNR mean (dB) for simulation.
        SimNoiseStd = 1;    % Default SNR standard deviation (dB) for simulation.
        SimDuration = 60;   % Default simulation duration (seconds).

        % --- Common Properties ---
        Fs = 1000;          % Sampling frequency (Hz). Assumed default, updated on data load.
        SelectedChannels = 1:16; % Indices of EEG channels to use.

        % --- Internal State ---
        rawTimeSeries % Cell array: {stimulus_idx/mean_idx, subject_idx/std_idx} -> data(sample, window, channel)
        currentExamTimeSeries % Temp storage for single loaded/generated exam data.
        isBulkLoaded = false % Flag indicating if bulk data is currently loaded.
        loadedStimulusIndices = [] % Indices of stimuli loaded in bulk mode.
        loadedSubjectIndices = []  % Indices of subjects loaded in bulk mode.
        simNoiseMeans = [] % Noise means used in bulk simulation.
        simNoiseStds = []  % Noise stds used in bulk simulation.
    end

    methods
        function obj = DataManager(mode, dataSourcePath, savePath)
            % Constructor for DataManager.
            % Args:
            %   mode: 'sim' or 'exp'.
            %   dataSourcePath: Path to experimental data or simulation settings.
            %   savePath: Path for saving results. Defaults to current directory.
            arguments
                mode {mustBeMember(mode, {'sim', 'exp'})}
                dataSourcePath string = ''
                savePath string = pwd % Default to current directory
            end
            obj.Mode = mode;
            obj.DataSourcePath = dataSourcePath;
            obj.SavePath = savePath;

            % Auto-generate noise frequencies for simulation if not provided
            if strcmp(obj.Mode, 'sim') && isempty(obj.SimNoiseFrequencies) && ~isempty(obj.SimSignalFrequencies)
                try
                    nNoise = numel(obj.SimSignalFrequencies);
                    poolSize = 200; % Size of the pool to draw from
                    freqStart = 300; % Starting frequency for noise pool
                    if poolSize >= nNoise
                       obj.SimNoiseFrequencies = freqStart + randperm(poolSize, nNoise);
                       fprintf('DataManager: Generated %d random noise frequencies (from %d Hz pool).\n', nNoise, freqStart);
                    else
                       warning('DataManager: Pool size (%d) is smaller than required noise frequencies (%d). Using replacement.', poolSize, nNoise);
                       obj.SimNoiseFrequencies = freqStart + randi(poolSize, 1, nNoise) -1; % Use replacement
                    end
                catch ME
                     warning('DataManager: Failed to auto-generate noise frequencies. %s',ME.message);
                end
            end
        end

        % --- Load Experimental Data Methods ---
        function obj = loadSingleExperiment(obj, subjectIdx, stimulusIdx)
            % Loads data for a single subject and stimulus.
            % Updates obj.currentExamTimeSeries and obj.Fs.
            if ~strcmp(obj.Mode, 'exp')
                error('DataManager: Cannot load experiment in simulation mode.');
            end
            if isempty(obj.DataSourcePath)
                error('DataManager: Experimental data source path is not set.');
            end
            if subjectIdx < 1 || subjectIdx > numel(obj.ExpSubjects) || stimulusIdx < 1 || stimulusIdx > numel(obj.ExpStimuli)
                error('DataManager: Invalid subject or stimulus index.');
            end

            subjectID = obj.ExpSubjects{subjectIdx};
            stimulusName = obj.ExpStimuli{stimulusIdx};
            fileName = [subjectID, stimulusName, '.mat'];
            filePath = fullfile(obj.DataSourcePath, fileName);

            try
                % Assuming .mat file contains 'x' (data) and 'Fs' (sampling rate)
                data = load(filePath, 'x', 'Fs');
                if ~isfield(data,'x') || ~isfield(data,'Fs')
                    error('DataManager: Loaded file %s missing required variables "x" or "Fs".', fileName);
                end
                obj.currentExamTimeSeries = data.x;
                obj.Fs = data.Fs; % Update Fs based on loaded data
                fprintf('DataManager: Loaded data for Subject %d (%s), Stimulus %d (%s). Fs = %d Hz.\n', ...
                        subjectIdx, subjectID, stimulusIdx, stimulusName, obj.Fs);
                obj.isBulkLoaded = false;
                obj.rawTimeSeries = {}; % Clear bulk data
            catch ME
                warning('DataManager: Failed to load file: %s\n  Error: %s', filePath, ME.message);
                obj.currentExamTimeSeries = []; % Indicate failure
            end
        end

        function obj = loadBulkExperiments(obj, subjectIndices, stimulusIndices)
            % Loads data for multiple subjects and stimuli specified by indices.
            % Stores data in obj.rawTimeSeries.
            if ~strcmp(obj.Mode, 'exp')
                error('DataManager: Cannot load experiments in simulation mode.');
            end
            if isempty(obj.DataSourcePath)
                error('DataManager: Experimental data source path is not set.');
            end
            if any(subjectIndices < 1) || any(subjectIndices > numel(obj.ExpSubjects)) || ...
               any(stimulusIndices < 1) || any(stimulusIndices > numel(obj.ExpStimuli))
                error('DataManager: Invalid subject or stimulus indices provided for bulk load.');
            end

            maxSubjIdx = max(subjectIndices);
            maxStimIdx = max(stimulusIndices);
            obj.rawTimeSeries = cell(maxStimIdx, maxSubjIdx); % Use max index for sizing

            fprintf('DataManager: Loading bulk experimental data...\n');
            loadedFs = []; % To check Fs consistency

            for stimIdx = stimulusIndices(:)' % Ensure row vector
                for subjIdx = subjectIndices(:)' % Ensure row vector
                    subjectID = obj.ExpSubjects{subjIdx};
                    stimulusName = obj.ExpStimuli{stimIdx};
                    fileName = [subjectID, stimulusName, '.mat'];
                    filePath = fullfile(obj.DataSourcePath, fileName);

                    try
                        data = load(filePath, 'x', 'Fs');
                         if ~isfield(data,'x') || ~isfield(data,'Fs')
                            error('DataManager: File %s missing "x" or "Fs".', fileName);
                        end
                        obj.rawTimeSeries{stimIdx, subjIdx} = data.x;

                        % Check Fs consistency
                        if isempty(loadedFs)
                            loadedFs = data.Fs;
                            obj.Fs = loadedFs; % Set Fs from first loaded file
                        elseif data.Fs ~= loadedFs
                             warning('DataManager: Inconsistent sampling rate found (File: %s, Fs=%d, Expected=%d). Using Fs=%d.', ...
                                     fileName, data.Fs, loadedFs, loadedFs);
                        end

                        fprintf('  Loaded: Subject %d (%s), Stimulus %d (%s)\n', subjIdx, subjectID, stimIdx, stimulusName);
                    catch ME
                        warning('  Failed to load: Subject %d (%s), Stimulus %d (%s)\n  Error: %s', subjIdx, subjectID, stimIdx, stimulusName, ME.message);
                        obj.rawTimeSeries{stimIdx, subjIdx} = []; % Mark as failed/empty
                    end
                end
            end
            obj.isBulkLoaded = true;
            obj.currentExamTimeSeries = []; % Clear single exam data
            obj.loadedStimulusIndices = stimulusIndices;
            obj.loadedSubjectIndices = subjectIndices;
            fprintf('DataManager: Bulk loading complete. Fs = %d Hz.\n', obj.Fs);
        end

        % --- Generate Simulated Data Methods ---
        function obj = generateSingleSimulation(obj, p)
            % Generates a single simulated dataset based on parameters.
            % Stores result in obj.currentExamTimeSeries.
             arguments
                obj
                p.duration double = obj.SimDuration
                p.noiseMean double = obj.SimNoiseMean
                p.noiseStd double = obj.SimNoiseStd
                p.signalFreqs double = obj.SimSignalFrequencies
                p.numChannels double = numel(obj.SelectedChannels)
                p.fs double = obj.Fs
            end
            if ~strcmp(obj.Mode, 'sim')
                error('DataManager: Cannot generate simulation in experimental mode.');
            end

            % --- Generate Base Signal ---
            totalSamplesPerChannel = round(p.fs * p.duration); % Ensure integer samples
            if totalSamplesPerChannel <= 0; error('DataManager: Duration results in non-positive sample count.'); end
            t = (0:totalSamplesPerChannel-1) / p.fs;

            signalBase = zeros(1, totalSamplesPerChannel);
            if ~isempty(p.signalFreqs)
                for freq = p.signalFreqs(:)' % Ensure row vector
                    signalBase = signalBase + sin(2*pi*freq*t);
                end
            else
                 fprintf('DataManager: No signal frequencies provided for simulation.\n');
            end


            % --- Apply Noise per Epoch (window) and Channel ---
            % Assume 1 epoch = 1 second = p.fs samples
            numEpochs = floor(totalSamplesPerChannel / p.fs);
            if numEpochs == 0; error('DataManager: Duration too short for any full epochs (1s).'); end
             samplesToUse = numEpochs * p.fs; % Use only full epochs

            simSignals = zeros(p.fs, numEpochs, p.numChannels); % samples_per_epoch x num_epochs x num_channels
            snrFun = @() p.noiseMean + p.noiseStd * randn();

            for chan = 1:p.numChannels
                % Reshape the base signal corresponding to full epochs
                 signalForEpochs = reshape(signalBase(1:samplesToUse), [p.fs, numEpochs]);

                for epoch = 1:numEpochs
                    epochSignal = signalForEpochs(:, epoch);
                    snr_dB = snrFun();
                    % AWGN adds noise based on signal power and desired SNR
                    noisySignal = awgn(epochSignal, snr_dB, 'measured', 'db');
                    % Scale to microvolts (typical EEG range)
                    simSignals(:, epoch, chan) = (10^-6) * noisySignal;
                end
            end
            obj.currentExamTimeSeries = simSignals;
            obj.isBulkLoaded = false;
            obj.rawTimeSeries = {};
             fprintf('DataManager: Generated single simulation (%d epochs, %d chan, Target SNR mean=%.1f dB).\n',...
                     numEpochs, p.numChannels, p.noiseMean);
         end

         function obj = generateBulkSimulations(obj, groupNoiseMeans, groupNoiseStds, p)
            % Generates multiple simulations varying noise parameters.
            % Stores results in obj.rawTimeSeries.
             arguments
                 obj
                 groupNoiseMeans double
                 groupNoiseStds double
                 p.duration double = obj.SimDuration
                 p.signalFreqs double = obj.SimSignalFrequencies
                 p.numChannels double = numel(obj.SelectedChannels)
                 p.fs double = obj.Fs
             end
              if ~strcmp(obj.Mode, 'sim')
                error('DataManager: Cannot generate simulations in experimental mode.');
              end

            numMeans = numel(groupNoiseMeans);
            numStds = numel(groupNoiseStds);
            obj.rawTimeSeries = cell(numMeans, numStds);

            fprintf('DataManager: Generating bulk simulation data (%d Means x %d Stds)...\n', numMeans, numStds);
            totalSims = numMeans * numStds;
            simCount = 0;

            for i = 1:numMeans
                meanVal = groupNoiseMeans(i);
                for j = 1:numStds
                    stdVal = groupNoiseStds(j);
                    simCount = simCount + 1;
                    % Generate a single simulation with these noise parameters
                    obj = obj.generateSingleSimulation(...
                        'duration', p.duration, ...
                        'noiseMean', meanVal, ...
                        'noiseStd', stdVal, ...
                        'signalFreqs', p.signalFreqs, ...
                        'numChannels', p.numChannels, ...
                        'fs', p.fs);
                    % Store the generated data in the bulk cell array
                    obj.rawTimeSeries{i, j} = obj.currentExamTimeSeries;
                    fprintf('  Generated Sim %d/%d: Mean=%.1f, Std=%.1f\n', simCount, totalSims, meanVal, stdVal);
                end
            end
            obj.isBulkLoaded = true;
            obj.currentExamTimeSeries = []; % Clear single exam temp storage
            obj.simNoiseMeans = groupNoiseMeans;
            obj.simNoiseStds = groupNoiseStds;
            fprintf('DataManager: Bulk simulation generation complete.\n');
         end

         % --- Utility and Data Access Methods ---
         function data = getTimeSeries(obj, idx1, idx2)
             % Returns data for a specific exam from the bulk-loaded set.
             % idx1: stimulus index (exp) or noise mean index (sim)
             % idx2: subject index (exp) or noise std index (sim)
             if ~obj.isBulkLoaded
                 warning('DataManager: Bulk data not loaded. Returning currentExamTimeSeries if available.');
                 data = obj.currentExamTimeSeries;
                 return;
             end
             if idx1 >= 1 && idx1 <= size(obj.rawTimeSeries, 1) && idx2 >= 1 && idx2 <= size(obj.rawTimeSeries, 2)
                 data = obj.rawTimeSeries{idx1, idx2};
             else
                 warning('DataManager: Indices (%d, %d) out of bounds for loaded bulk data [%d x %d].', ...
                         idx1, idx2, size(obj.rawTimeSeries, 1), size(obj.rawTimeSeries, 2));
                 data = [];
             end
         end

         function bulkData = getBulkTimeSeries(obj)
             % Returns the entire cell array of bulk-loaded raw data.
             if ~obj.isBulkLoaded
                 warning('DataManager: Bulk data not loaded. Returning empty cell.');
                 bulkData = {};
             else
                bulkData = obj.rawTimeSeries;
             end
         end

         function [stimIndices, subjIndices] = getLoadedIndices(obj)
             % Returns the indices used during the last bulk load/generation.
             if ~obj.isBulkLoaded
                  stimIndices = []; subjIndices = [];
                  warning('DataManager: No bulk data loaded, returning empty indices.');
             elseif strcmp(obj.Mode,'exp')
                 stimIndices = obj.loadedStimulusIndices;
                 subjIndices = obj.loadedSubjectIndices;
             else % sim
                 stimIndices = 1:numel(obj.simNoiseMeans); % Represent mean index
                 subjIndices = 1:numel(obj.simNoiseStds); % Represent std index
             end
         end

         function nBins = getNumBins(obj)
            % Calculate number of FFT bins (positive frequencies).
            nBins = floor(obj.Fs / 2) + 1;
         end

         function nfft = getNFFT(obj)
             % Number of FFT points (often equal to Fs for 1-sec resolution).
             nfft = obj.Fs;
         end

         function Mmax = getSuggestedMMax(obj, stimulusIdx)
             % Returns the suggested Mmax for a given experimental stimulus index.
             if strcmp(obj.Mode, 'exp') && stimulusIdx >= 1 && stimulusIdx <= numel(obj.ExpSuggestedMMax)
                 Mmax = obj.ExpSuggestedMMax(stimulusIdx);
             else
                 Mmax = []; % Not applicable or index out of range
             end
         end

    end % methods
end % classdef