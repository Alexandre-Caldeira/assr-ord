classdef PreProcessor < handle
    % Applies preprocessing steps to raw time-domain EEG signals.

    properties (SetAccess = protected)
        zanoteliGain = 200;
        zanoteliCutThreshold
        filterFcLower = 70;
        filterFcUpper = 499;
        filterOrder = 8;
        fs = 1000;
    end

    methods
        function obj = PreProcessor(options)
            arguments
                options.zanoteliGain (1,1) {mustBeNumeric, mustBePositive} = 200
                options.filterFcLower (1,1) {mustBeNumeric, mustBePositive} = 70
                options.filterFcUpper (1,1) {mustBeNumeric, mustBePositive} = 499
                options.filterOrder (1,1) {mustBeInteger, mustBePositive} = 8
                options.fs (1,1) {mustBeNumeric, mustBePositive} = 1000
            end

            obj.zanoteliGain = options.zanoteliGain;
            obj.filterFcLower = options.filterFcLower;
            obj.filterFcUpper = options.filterFcUpper;
            obj.filterOrder = options.filterOrder;
            obj.fs = options.fs;
            obj.zanoteliCutThreshold = 0.1 / obj.zanoteliGain;

            nyquist = obj.fs / 2;
            if obj.filterFcUpper >= nyquist
                obj.filterFcUpper = nyquist - 1;
                 warning('Filter upper cutoff frequency adjusted to %.1f Hz (Nyquist).', obj.filterFcUpper);
            end
            if obj.filterFcLower >= obj.filterFcUpper
                warning('Filter lower cutoff >= upper cutoff. Disabling filter.', obj.filterFcLower, obj.filterFcUpper);
                obj.filterFcLower = -1;
            end
            % Reduced verbosity
            % fprintf('[%s] PreProcessor initialized. Zanoteli Threshold: %.4e. Filter: [%.1f, %.1f] Hz, Order: %d.\n', ...
            %        datetime, obj.zanoteliCutThreshold, obj.filterFcLower, obj.filterFcUpper, obj.filterOrder);
        end

        function processedSignals = applyZanoteliPreprocessing(obj, rawSignals, maxWindows)
            arguments
                obj
                rawSignals (:,:,:) {mustBeNumeric}
                maxWindows (1,1) {mustBeInteger, mustBePositive} = size(rawSignals, 2)
            end

            nSamples = size(rawSignals, 1);
            nChannels = size(rawSignals, 3);
            if nSamples ~= obj.fs
                warning('Input signal Fs (%d samples/window) does not match Preprocessor Fs (%d).', nSamples, obj.fs);
            end

            processedSignals = zeros(nSamples, maxWindows, nChannels);
            finalWindowCount = maxWindows;

            % Reduced verbosity
            % fprintf('[%s] Applying Zanoteli Preprocessing...\n', datetime);

            for ch = 1:nChannels
                x_ch = squeeze(rawSignals(:, :, ch));
                x_ch = x_ch - mean(x_ch, 1);
                if size(x_ch, 2) > 2, x_ch = x_ch(:, 3:end); end
                if isempty(x_ch), finalWindowCount = 0; continue; end

                maxAbsPerWindow = max(abs(x_ch), [], 1);
                artifactIndices = find(maxAbsPerWindow > obj.zanoteliCutThreshold);
                validIndices = setdiff(1:size(x_ch, 2), artifactIndices);
                x_clean = x_ch(:, validIndices);

                numValidWindows = size(x_clean, 2);
                windowsToKeep = min(numValidWindows, maxWindows);

                if windowsToKeep > 0
                   processedSignals(:, 1:windowsToKeep, ch) = x_clean(:, 1:windowsToKeep);
                end
                finalWindowCount = min(finalWindowCount, windowsToKeep);
            end

            if finalWindowCount < maxWindows || finalWindowCount == 0
                if finalWindowCount > 0
                    processedSignals = processedSignals(:, 1:finalWindowCount, :);
                    % Reduced verbosity
                    % fprintf('[%s] Zanoteli Preprocessing complete. Final signal size: %d samples x %d windows x %d channels.\n', ...
                    %        datetime, size(processedSignals,1), size(processedSignals,2), size(processedSignals,3));
                else
                    processedSignals = [];
                     % Reduced verbosity
                     % fprintf('[%s] Zanoteli Preprocessing complete. No valid windows remained.\n', datetime);
                end
            % else % Reduced verbosity
                 % fprintf('[%s] Zanoteli Preprocessing complete. Kept max %d windows.\n', datetime, maxWindows);
            end
        end

        function filteredSignals = applyAntunesFiltering(obj, processedSignals)
            arguments
                obj
                processedSignals (:,:,:) {mustBeNumeric}
            end

            if obj.filterFcLower < 0
                 warning('Filtering disabled due to invalid cutoff frequencies.');
                 filteredSignals = processedSignals; return;
            end
            if size(processedSignals, 2) < obj.filterOrder * 3
                 warning('Skipping filtering: Signal window count (%d) too short for filter order (%d).', size(processedSignals, 2), obj.filterOrder);
                 filteredSignals = processedSignals; return;
            end

            % Reduced verbosity
            % fprintf('[%s] Applying Antunes Filtering ([%.1f, %.1f] Hz, Order %d)...\n', ...
            %        datetime, obj.filterFcLower, obj.filterFcUpper, obj.filterOrder);

            filteredSignals = zeros(size(processedSignals));
            nChannels = size(processedSignals, 3);
            nyquist = obj.fs / 2;
            Wn = [obj.filterFcLower, obj.filterFcUpper] / nyquist;

            try, [b, a] = butter(obj.filterOrder, Wn, 'bandpass');
            catch ME, warning('Could not design Butterworth filter: %s. Skipping filtering.', ME.message); filteredSignals = processedSignals; return; end

            for ch = 1:nChannels
                x_ch = squeeze(processedSignals(:, :, ch));
                filtered_x_ch = zeros(size(x_ch));
                for win = 1:size(x_ch, 2)
                   try, filtered_x_ch(:, win) = filtfilt(b, a, x_ch(:, win));
                   catch ME_filt, warning('Filtfilt error on chan %d, win %d: %s. Skipping window.', ch, win, ME_filt.message); filtered_x_ch(:, win) = x_ch(:, win); end % Handle potential filtfilt errors
                end
                filteredSignals(:, :, ch) = filtered_x_ch;
            end
             % Reduced verbosity
             % fprintf('[%s] Antunes Filtering complete.\n', datetime);
        end

        function processedSignals = processSignals(obj, rawSignals, maxWindows, useFilter)
             arguments
                obj
                rawSignals (:,:,:) {mustBeNumeric}
                maxWindows (1,1) {mustBeInteger, mustBePositive} = size(rawSignals, 2)
                useFilter (1,1) logical = true % Control filtering step
            end
            processedSignals = obj.applyZanoteliPreprocessing(rawSignals, maxWindows);
            if useFilter && ~isempty(processedSignals)
                processedSignals = obj.applyAntunesFiltering(processedSignals);
            end
        end
    end
end