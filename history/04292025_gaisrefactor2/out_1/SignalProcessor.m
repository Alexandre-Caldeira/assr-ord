% --- FILE START: SignalProcessor.m ---
classdef SignalProcessor
    % Handles time-domain preprocessing (artifact rejection, filtering)
    % and FFT calculation for EEG data.

    properties
        % --- Zanoteli Preprocessing Parameters ---
        ZanoteliEnable (1,1) logical = true
        ZanoteliGain (1,1) double {mustBePositive} = 200 % Keep validation for independent params
        ZanoteliThreshold (1,1) double % Calculated

        % --- Filtering Parameters ---
        FilterEnable (1,1) logical = true
        FilterFcLower (1,1) double {mustBePositive} = 70 % Keep validation
        FilterFcUpper (1,1) double       % Calculated
        FilterOrder (1,1) double {mustBeInteger, mustBePositive} = 8 % Keep validation

        % --- Output Storage ---
        processedTimeSeries % Cell array {stim/mean_idx, subj/std_idx} -> data(sample, window, channel)
        frequencyDomainData % Cell array {stim/mean_idx, subj/std_idx} -> FFT data(bin, window, channel)

        % --- Internal State ---
        % Remove validation from properties derived from DataManager
        Fs (1,1) double
        NFFT (1,1) double
        NumBins (1,1) double
        NumChannels (1,1) double
    end

    methods
        function obj = SignalProcessor(dataManager)
            % Constructor: Initializes processor with info from DataManager.
            % All properties derived from dataManager are calculated AND validated here.
            arguments
                dataManager DataManager
            end

            % --- 1. Get DataManager Properties ---
            tempFs = dataManager.Fs;
            tempNFFT = dataManager.getNFFT(); % Should be equal to Fs usually
            tempNumBins = dataManager.getNumBins();
            tempNumChannels = numel(dataManager.SelectedChannels);

            % --- 2. Validate Retrieved Properties ---
            if ~isscalar(tempFs) || tempFs <= 0
                error('SignalProcessor: Fs value (%g) received from DataManager must be a positive scalar.', tempFs);
            end
             if ~isscalar(tempNFFT) || tempNFFT <= 0 || floor(tempNFFT) ~= tempNFFT
                error('SignalProcessor: NFFT value (%g) derived from DataManager must be a positive integer scalar.', tempNFFT);
            end
             if ~isscalar(tempNumBins) || tempNumBins <= 0 || floor(tempNumBins) ~= tempNumBins
                error('SignalProcessor: NumBins value (%g) derived from DataManager must be a positive integer scalar.', tempNumBins);
            end
             if ~isscalar(tempNumChannels) || tempNumChannels <= 0 || floor(tempNumChannels) ~= tempNumChannels
                 error('SignalProcessor: NumChannels value (%g) derived from DataManager must be a positive integer scalar.', tempNumChannels);
            end

            % --- 3. Assign Validated Properties to Object ---
            obj.Fs = tempFs;
            obj.NFFT = tempNFFT;
            obj.NumBins = tempNumBins;
            obj.NumChannels = tempNumChannels;

            % --- 4. Calculate and Validate Dependent Properties ---
            obj.ZanoteliThreshold = 0.1 / obj.ZanoteliGain; % Threshold in Volts
            obj.FilterFcUpper = obj.Fs / 2 - 1; % Nyquist minus a bit

            if obj.FilterFcUpper <= 0
                error('SignalProcessor: Calculated FilterFcUpper (%g) must be positive. Check Fs value (%g).', obj.FilterFcUpper, obj.Fs);
            end
             if obj.FilterEnable && obj.FilterFcLower >= obj.FilterFcUpper
                warning('SignalProcessor: Lower filter cutoff (%g Hz) is not less than upper cutoff (%g Hz). Adjust parameters or disable filter.', obj.FilterFcLower, obj.FilterFcUpper);
             end

            fprintf('SignalProcessor: Initialized. Fs=%d, NFFT=%d, NumBins=%d, NumChannels=%d.\n', obj.Fs, obj.NFFT, obj.NumBins, obj.NumChannels);
        end % Constructor

        % --- applyPreprocessing Method ---
        function processedSignal = applyPreprocessing(obj, rawSignal, Mmax_limit)
             % Applies configured preprocessing steps (Zanoteli, Filtering) to a single exam's data.
             % Args:
             %   rawSignal: Raw time series data (samples x windows x channels).
             %   Mmax_limit: Maximum number of windows to retain after artifact rejection.
             % Returns:
             %   processedSignal: Preprocessed time series data, potentially shorter than Mmax_limit. Returns empty if processing fails.

             processedSignal = []; % Initialize output
             if isempty(rawSignal)
                 fprintf('applyPreprocessing: Input rawSignal is empty.\n');
                 return;
             end

             currentSignal = rawSignal;
             [nSamples, nWindowsRaw, nChanRaw] = size(currentSignal);
             if nChanRaw == 0; fprintf('applyPreprocessing: Input rawSignal has 0 channels.\n'); return; end
             if nSamples == 0; fprintf('applyPreprocessing: Input rawSignal has 0 samples.\n'); return; end
             if nWindowsRaw == 0; fprintf('applyPreprocessing: Input rawSignal has 0 windows.\n'); return; end


             % Use the actual number of channels from the input data for processing this specific signal
             numChanToProcess = nChanRaw;
             if nChanRaw ~= obj.NumChannels
                warning('applyPreprocessing: Channel count in data (%d) differs from processor config (%d). Processing %d channels.', nChanRaw, obj.NumChannels, nChanRaw);
             end

             % --- 1. Zanoteli Artifact Removal (if enabled) ---
             if obj.ZanoteliEnable
                 % Use Mmax_limit which should be Inf if not specified or calculable
                 effectiveMmax = min(Mmax_limit, nWindowsRaw); % Don't preallocate beyond original windows unless necessary
                 if isinf(effectiveMmax); effectiveMmax = nWindowsRaw; end % Handle Inf case

                 fprintf('  Applying Zanoteli preprocessing (Threshold=%.2e V, Mmax Limit=%d)...\n', obj.ZanoteliThreshold, effectiveMmax);
                 processedZanoteli = nan(nSamples, effectiveMmax, numChanToProcess); % Preallocate size based on Mmax limit or original size
                 minWindowsAcrossChannels = effectiveMmax;

                 for chan = 1:numChanToProcess
                     x_chan = squeeze(currentSignal(:, :, chan));

                     % Remove DC per window (epoch)
                     try % Add try-catch for robustness
                        x_chan_mean = mean(x_chan, 1);
                        x_chan = x_chan - x_chan_mean;
                     catch ME_dc
                         warning('Zanoteli: Error removing DC for channel %d: %s. Skipping channel.', chan, ME_dc.message);
                         minWindowsAcrossChannels = 0; continue;
                     end

                     % Skip first 2 seconds (windows)
                     if size(x_chan, 2) > 2
                         x_chan(:, 1:2) = [];
                     else
                         warning('Zanoteli: Signal too short (<3 windows) to remove first 2s for channel %d. Skipping channel.', chan);
                         minWindowsAcrossChannels = 0; continue;
                     end

                     % Amplitude thresholding
                     try
                        Vmax = max(abs(x_chan), [], 1);
                        valid_windows_mask = Vmax <= obj.ZanoteliThreshold;
                        x_clean = x_chan(:, valid_windows_mask);
                     catch ME_thresh
                          warning('Zanoteli: Error during thresholding for channel %d: %s. Skipping channel.', chan, ME_thresh.message);
                          minWindowsAcrossChannels = 0; continue;
                     end

                     % Limit to Mmax_limit and store valid portion
                     actualValidWindows = size(x_clean, 2);
                     windowsToKeep = min(actualValidWindows, effectiveMmax); % Limit by original effectiveMmax

                     if windowsToKeep > 0
                         % Only assign the kept windows to the correctly sized preallocated array
                         processedZanoteli(:, 1:windowsToKeep, chan) = x_clean(:, 1:windowsToKeep);
                         minWindowsAcrossChannels = min(minWindowsAcrossChannels, windowsToKeep);
                     else
                         warning('  Zanoteli: No windows below threshold for channel %d after removing first 2.', chan);
                         minWindowsAcrossChannels = 0; % No valid windows found for this channel AT ALL
                         % Break inner loop? If one channel has 0 valid, maybe whole dataset invalid? Depends on requirement.
                         % Assuming we continue and trim later.
                     end
                 end % End channel loop

                 % Trim all channels to the minimum common length
                  if minWindowsAcrossChannels > 0
                     % Trim the NaN preallocated array to the actual min common length
                     currentSignal = processedZanoteli(:, 1:minWindowsAcrossChannels, :);
                      % Check if trimming resulted in empty signal (should only happen if minWindows=0)
                      if isempty(currentSignal) || size(currentSignal,2)==0
                          warning('  Zanoteli: Trimming resulted in 0 common valid windows. Returning empty.');
                          processedSignal = []; return;
                      end
                     fprintf('  Zanoteli complete. Resulting windows: %d.\n', minWindowsAcrossChannels);
                 else
                     warning('  Zanoteli: Preprocessing resulted in 0 common valid windows across channels. Returning empty.');
                     processedSignal = []; return;
                 end

             else % Zanoteli disabled - apply Mmax_limit if it's finite
                 if ~isinf(Mmax_limit)
                     nWindowsCurrent = size(currentSignal, 2);
                     windowsToKeep = min(nWindowsCurrent, Mmax_limit);
                     if windowsToKeep < nWindowsCurrent
                         currentSignal = currentSignal(:, 1:windowsToKeep, :);
                         fprintf('  Preprocessing: Signal trimmed to Mmax limit %d.\n', windowsToKeep);
                     end
                 end
             end % End ZanoteliEnable check

             % --- Signal Empty Check ---
             if isempty(currentSignal) || size(currentSignal,2) == 0
                  processedSignal = [];
                  fprintf('applyPreprocessing: Signal is empty after initial processing/trimming steps.\n');
                  return;
             end

             % --- 2. Filtering ---
             if obj.FilterEnable
                 [nSamplesFilt, nWindowsFilt, nChanFilt] = size(currentSignal);
                 fprintf('  Applying Butterworth filter (Order=%d, Fc=[%.1f, %.1f] Hz)...\n', obj.FilterOrder, obj.FilterFcLower, obj.FilterFcUpper);

                 filtfiltMinLength = obj.FilterOrder * 3 + 1; % filtfilt needs length > 3*order
                 if nSamplesFilt < filtfiltMinLength
                     warning('  Filtering: Signal length per window (%d samples) is too short for stable filtfilt (Order %d requires >=%d). Skipping filter.', nSamplesFilt, obj.FilterOrder, filtfiltMinLength);
                     processedSignal = currentSignal; return;
                 end

                 try
                    [b, a] = butter(obj.FilterOrder, [obj.FilterFcLower, obj.FilterFcUpper] / (obj.Fs / 2), 'bandpass');
                 catch ME_butter
                      warning('Filtering: Could not design Butterworth filter. Error: %s. Skipping filter.', ME_butter.message);
                      processedSignal = currentSignal; return;
                 end

                 filteredSignal = zeros(size(currentSignal)); % Preallocate
                 for chan = 1:nChanFilt
                     signal_chan_win = squeeze(currentSignal(:, :, chan));
                     try
                         filtered_chan_win = filtfilt(b, a, signal_chan_win); % Apply column-wise
                         filteredSignal(:, :, chan) = filtered_chan_win;
                     catch ME_filt
                          warning('Filtering: Error during filtfilt on channel %d: %s. Keeping original signal for this channel.', chan, ME_filt.message);
                          filteredSignal(:, :, chan) = signal_chan_win; % Keep original
                     end
                 end
                 currentSignal = filteredSignal;
                 fprintf('  Filtering complete.\n');
             end % End FilterEnable check

             processedSignal = currentSignal;
             fprintf('applyPreprocessing: Finished. Output size: %d x %d x %d\n', size(processedSignal));
        end % applyPreprocessing


         % --- applyFFT Method ---
         function fftSignal = applyFFT(obj, processedSignal)
            % Applies FFT to preprocessed time-domain data (samples x windows x channels).
            % Returns FFT data (bins x windows x channels).

            fftSignal = []; % Initialize output
            if isempty(processedSignal)
                 fprintf('applyFFT: Input processedSignal is empty.\n');
                return;
            end

            [nSamples, nWindows, nChan] = size(processedSignal);
            if nChan == 0 || nWindows == 0 || nSamples == 0
                fprintf('applyFFT: Input processedSignal has zero dimension (%d x %d x %d).\n', nSamples, nWindows, nChan);
                return;
            end

             numChanToProcess = nChan;
             if nChan ~= obj.NumChannels
                warning('applyFFT: Channel count in data (%d) differs from processor config (%d). Processing %d channels.', nChan, obj.NumChannels, nChan);
             end

            fftSignal = zeros(obj.NumBins, nWindows, numChanToProcess); % Preallocate

            fprintf('  Applying FFT (NFFT=%d, NumBins=%d)...\n', obj.NFFT, obj.NumBins);
            for chan = 1:numChanToProcess
                for epoch = 1:nWindows
                    epochData = squeeze(processedSignal(:, epoch, chan));
                    if any(isnan(epochData))
                        warning('  FFT: NaNs found in data for epoch %d, channel %d. FFT result will be NaN.', epoch, chan);
                        temp = nan(obj.NFFT, 1);
                    else
                         try
                            temp = fft(epochData, obj.NFFT);
                         catch ME_fft
                             warning('  FFT: Error during fft for epoch %d, channel %d: %s. Result set to NaN.', epoch, chan, ME_fft.message);
                             temp = nan(obj.NFFT, 1);
                         end
                    end
                    % Check if temp is correct size before indexing
                    if numel(temp) >= obj.NumBins
                        fftSignal(:, epoch, chan) = temp(1:obj.NumBins);
                    else
                         warning('  FFT: Result size (%d) smaller than NumBins (%d) for epoch %d, channel %d. Padding with NaN.', numel(temp), obj.NumBins, epoch, chan);
                         fftSignal(:, epoch, chan) = NaN; % Mark as invalid
                    end
                end
            end
             fprintf('  FFT complete.\n');
         end % applyFFT

         % --- processBulkData Method ---
         function obj = processBulkData(obj, dataManager)
            % Processes bulk data stored in DataManager through preprocessing and FFT.
            % Stores results in obj.processedTimeSeries and obj.frequencyDomainData.
            arguments
                 obj SignalProcessor
                 dataManager DataManager
            end

            if ~dataManager.isBulkLoaded
                error('SignalProcessor: Bulk data not loaded in DataManager. Load data first.');
            end

            bulkRawData = dataManager.getBulkTimeSeries();
            if isempty(bulkRawData)
                error('SignalProcessor: Bulk data retrieved from DataManager is empty.');
            end

            [nDim1, nDim2] = size(bulkRawData);
            obj.processedTimeSeries = cell(nDim1, nDim2);
            obj.frequencyDomainData = cell(nDim1, nDim2);

            fprintf('SignalProcessor: Starting bulk processing (%d x %d exams)...\n', nDim1, nDim2);
            totalExams = numel(bulkRawData);
            processedCount = 0;

            for idx1 = 1:nDim1 % Iterate stimulus/mean index
                 % Determine Mmax limit for this primary index (Stimulus or Mean Group)
                 suggestedMmax = Inf; % Default
                 max_len_found = 0; % Track max windows found for this index1 group
                 needs_mmax_calc = false;

                 if strcmp(dataManager.Mode,'exp')
                    mmax_val = dataManager.getSuggestedMMax(idx1);
                    if ~isempty(mmax_val)
                        suggestedMmax = mmax_val;
                    else
                        needs_mmax_calc = true; % Need to find max length across subjects
                        warning('SignalProcessor:processBulkData', 'No ExpSuggestedMMax for Stimulus Index %d.', idx1);
                    end
                 else % Simulation mode
                      needs_mmax_calc = true; % Always find max length for sims unless duration was strictly enforced previously
                 end

                 % If needed, find the maximum number of windows across the secondary index for this primary index
                 if needs_mmax_calc
                      max_windows_group = 0;
                      for idx2_inner = 1:nDim2
                          if idx1 <= size(bulkRawData,1) && idx2_inner <= size(bulkRawData,2) && ~isempty(bulkRawData{idx1, idx2_inner})
                              max_windows_group = max(max_windows_group, size(bulkRawData{idx1, idx2_inner}, 2));
                          end
                      end
                      if max_windows_group > 0
                           suggestedMmax = max_windows_group;
                           fprintf('  Setting Mmax limit based on observed data for Index1=%d: %d windows.\n', idx1, suggestedMmax);
                      else
                           warning('SignalProcessor:processBulkData','No data found to determine Mmax limit for Index1=%d. Using Inf.', idx1);
                           suggestedMmax = Inf; % Fallback
                      end
                 end
                 % End Mmax determination for this group

                for idx2 = 1:nDim2 % Iterate subject/std index
                    processedCount = processedCount + 1;
                    fprintf('Processing Exam %d/%d (Index1: %d, Index2: %d)...\n', processedCount, totalExams, idx1, idx2);
                    rawData = []; % Ensure rawData is reset or defined
                    if idx1 <= size(bulkRawData,1) && idx2 <= size(bulkRawData,2)
                        rawData = bulkRawData{idx1, idx2};
                    end

                    if isempty(rawData)
                        fprintf('  Skipping: No raw data available for Index1=%d, Index2=%d.\n', idx1, idx2);
                         obj.processedTimeSeries{idx1, idx2} = [];
                         obj.frequencyDomainData{idx1, idx2} = [];
                        continue;
                    end

                    % --- 1. Preprocess ---
                    processedData = obj.applyPreprocessing(rawData, suggestedMmax);
                    obj.processedTimeSeries{idx1, idx2} = processedData;

                    % --- 2. FFT ---
                    fftData = obj.applyFFT(processedData); % Handles empty processedData
                    obj.frequencyDomainData{idx1, idx2} = fftData;

                     fprintf('  Finished Exam %d/%d.\n', processedCount, totalExams);

                end % End inner loop (idx2)
            end % End outer loop (idx1)
            fprintf('SignalProcessor: Bulk processing complete.\n');
         end % processBulkData

         % --- Access Processed Data ---
         function data = getProcessedTimeSeries(obj, idx1, idx2)
            % Access preprocessed time series for a specific exam.
             if idx1 >= 1 && idx1 <= size(obj.processedTimeSeries, 1) && idx2 >= 1 && idx2 <= size(obj.processedTimeSeries, 2)
                 data = obj.processedTimeSeries{idx1, idx2};
             else
                 warning('SignalProcessor:getProcessedTimeSeries', 'Indices (%d, %d) out of bounds for processed time series data [%d x %d].', idx1, idx2, size(obj.processedTimeSeries,1), size(obj.processedTimeSeries,2));
                 data = [];
             end
         end

         function data = getFrequencyDomainData(obj, idx1, idx2)
            % Access frequency domain (FFT) data for a specific exam.
             if idx1 >= 1 && idx1 <= size(obj.frequencyDomainData, 1) && idx2 >= 1 && idx2 <= size(obj.frequencyDomainData, 2)
                 data = obj.frequencyDomainData{idx1, idx2};
             else
                  warning('SignalProcessor:getFrequencyDomainData', 'Indices (%d, %d) out of bounds for frequency domain data [%d x %d].', idx1, idx2, size(obj.frequencyDomainData,1), size(obj.frequencyDomainData,2));
                 data = [];
             end
         end

         function bulkData = getBulkFrequencyDomainData(obj)
            % Returns the entire cell array of frequency domain data.
             bulkData = obj.frequencyDomainData;
         end

         function bulkData = getBulkProcessedTimeSeries(obj)
             % Returns the entire cell array of processed time series data.
             bulkData = obj.processedTimeSeries;
         end

    end % methods
end % classdef