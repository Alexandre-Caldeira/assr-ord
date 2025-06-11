classdef PreProcessor < handle
    % Applies preprocessing steps to raw time-domain EEG signals.
    % Includes artifact rejection (Zanoteli method) and bandpass filtering.

    properties (SetAccess = protected)
        % Zanoteli Artifact Rejection Parameters
        zanoteliGain = 200;
        zanoteliCutThreshold % Calculated based on gain

        % Butterworth Filter Parameters
        filterFcLower = 70;  % Lower cutoff frequency (Hz)
        filterFcUpper = 499; % Upper cutoff frequency (Hz), ensure below Nyquist
        filterOrder = 8;     % Butterworth filter order

        fs = 1000; % Sampling frequency - MUST match data
    end

    methods
        function obj = PreProcessor(options)
            % Constructor for PreProcessor
            % options: Name-Value pairs:
            %   'zanoteliGain': Gain for artifact rejection threshold
            %   'filterFcLower': Lower cutoff Hz
            %   'filterFcUpper': Upper cutoff Hz
            %   'filterOrder': Filter order
            %   'fs': Sampling frequency (MUST match data)

            arguments % Define options block for Name-Value pairs
                options.zanoteliGain (1,1) {mustBeNumeric, mustBePositive} = 200
                options.filterFcLower (1,1) {mustBeNumeric, mustBePositive} = 70
                options.filterFcUpper (1,1) {mustBeNumeric, mustBePositive} = 499
                options.filterOrder (1,1) {mustBeInteger, mustBePositive} = 8
                options.fs (1,1) {mustBeNumeric, mustBePositive} = 1000
            end

            % Assign properties directly from the parsed options
            obj.zanoteliGain = options.zanoteliGain;
            obj.filterFcLower = options.filterFcLower;
            obj.filterFcUpper = options.filterFcUpper;
            obj.filterOrder = options.filterOrder;
            obj.fs = options.fs;

             % Calculate Zanoteli threshold
             obj.zanoteliCutThreshold = 0.1 / obj.zanoteliGain;

             % Validate filter upper cutoff against Nyquist
             nyquist = obj.fs / 2;
             if obj.filterFcUpper >= nyquist
                 warning('Filter upper cutoff frequency (%.1f Hz) is >= Nyquist frequency (%.1f Hz). Adjusting to %.1f Hz.', ...
                         obj.filterFcUpper, nyquist, nyquist - 1);
                 obj.filterFcUpper = nyquist - 1;
             end
             if obj.filterFcLower >= obj.filterFcUpper
                  warning('Filter lower cutoff frequency (%.1f Hz) is >= upper cutoff (%.1f Hz). Disabling filter.', ...
                         obj.filterFcLower, obj.filterFcUpper);
                 obj.filterFcLower = -1; % Flag to disable filter
             end
             fprintf('[%s] PreProcessor initialized. Zanoteli Threshold: %.4e. Filter: [%.1f, %.1f] Hz, Order: %d.\n', ...
                     datetime, obj.zanoteliCutThreshold, obj.filterFcLower, obj.filterFcUpper, obj.filterOrder);
        end

                function processedSignals = applyZanoteliPreprocessing(obj, rawSignals, maxWindows)
            % Applies Zanoteli's artifact rejection and DC removal.
            % rawSignals: Input signal matrix (samples, windows, channels)
            % maxWindows: Maximum number of windows to keep after rejection
            % returns: Processed signal matrix (potentially fewer windows)

            arguments
                obj
                rawSignals (:,:,:) {mustBeNumeric}
                maxWindows (1,1) {mustBeInteger, mustBePositive} = size(rawSignals, 2) % Default: keep all possible
            end

            % --- Get input dimensions ---
            nSamples = size(rawSignals, 1); % Use ACTUAL samples from input
            nChannels = size(rawSignals, 3);
            % ---

            % --- Check against stored Fs (as a warning) ---
            if nSamples ~= obj.fs
                warning('Input signal Fs (%d samples/window) does not match Preprocessor Fs (%d). Filter calculations may be affected if different Fs intended.', ...
                        nSamples, obj.fs);
            end
            % ---

            % --- Preallocate based on INPUT samples ---
            processedSignals = zeros(nSamples, maxWindows, nChannels); % Use nSamples here
            % ---
            finalWindowCount = maxWindows; % Track actual windows kept

            fprintf('[%s] Applying Zanoteli Preprocessing...\n', datetime);

            for ch = 1:nChannels
                x_ch = squeeze(rawSignals(:, :, ch)); % Samples x Windows

                % 1. Remove DC offset per window (epoch)
                x_ch = x_ch - mean(x_ch, 1); % Subtract column means

                % 2. Exclude first 2 seconds (windows) if available
                if size(x_ch, 2) > 2
                    x_ch = x_ch(:, 3:end);
                % else % No warning needed if signal just short
                %    warning('Signal too short (< 3 windows) to remove first 2s for channel %d.', ch);
                end

                % Check if any windows remain after step 2
                if isempty(x_ch)
                    finalWindowCount = 0; % No windows left for this channel
                    continue; % Skip to next channel if no data remains
                end

                % 3. Artifact rejection based on amplitude
                maxAbsPerWindow = max(abs(x_ch), [], 1);
                artifactIndices = find(maxAbsPerWindow > obj.zanoteliCutThreshold);
                validIndices = setdiff(1:size(x_ch, 2), artifactIndices);

                x_clean = x_ch(:, validIndices);

                % 4. Limit number of windows
                numValidWindows = size(x_clean, 2);
                windowsToKeep = min(numValidWindows, maxWindows);

                if windowsToKeep < maxWindows && numValidWindows > 0 % Only report if some were kept but fewer than max
                    fprintf('  Channel %d: Kept %d windows after artifact rejection (target max: %d).\n', ...
                            ch, windowsToKeep, maxWindows);
                end

                % Store in preallocated matrix, handle case where fewer windows remain
                if windowsToKeep > 0 % Ensure there's something to copy
                   processedSignals(:, 1:windowsToKeep, ch) = x_clean(:, 1:windowsToKeep); % Assignment should now work
                end

                % Update the final window count based on the minimum across channels
                finalWindowCount = min(finalWindowCount, windowsToKeep);
            end

            % Trim the output matrix if fewer than maxWindows were kept consistently
            % across all channels, or if some channels had zero windows left.
            if finalWindowCount < maxWindows || finalWindowCount == 0
                if finalWindowCount > 0
                    processedSignals = processedSignals(:, 1:finalWindowCount, :);
                    fprintf('[%s] Zanoteli Preprocessing complete. Final signal size: %d samples x %d windows x %d channels.\n', ...
                            datetime, size(processedSignals,1), size(processedSignals,2), size(processedSignals,3));
                else
                    processedSignals = []; % Return empty if no windows survived
                     fprintf('[%s] Zanoteli Preprocessing complete. No valid windows remained after rejection.\n', datetime);
                end
            else
                 fprintf('[%s] Zanoteli Preprocessing complete. Kept max %d windows.\n', datetime, maxWindows);
            end

        end % End applyZanoteliPreprocessing

        function filteredSignals = applyAntunesFiltering(obj, processedSignals)
            % Applies Butterworth bandpass filter to processed signals.
            arguments
                obj
                processedSignals (:,:,:) {mustBeNumeric}
            end

            if obj.filterFcLower < 0 % Check if filtering is disabled
                 warning('Filtering disabled due to invalid cutoff frequencies.');
                 filteredSignals = processedSignals;
                 return
            end
            
            if size(processedSignals, 2) < obj.filterOrder * 3 % filtfilt requirement
                 warning('Skipping filtering: Signal window count (%d) too short for filter order (%d). Needs > %d.', ...
                         size(processedSignals, 2), obj.filterOrder, obj.filterOrder * 3);
                 filteredSignals = processedSignals;
                 return
            end

            fprintf('[%s] Applying Antunes Filtering ([%.1f, %.1f] Hz, Order %d)...\n', ...
                    datetime, obj.filterFcLower, obj.filterFcUpper, obj.filterOrder);

            filteredSignals = zeros(size(processedSignals));
            nChannels = size(processedSignals, 3);
            nyquist = obj.fs / 2;
            Wn = [obj.filterFcLower, obj.filterFcUpper] / nyquist;

            try
                [b, a] = butter(obj.filterOrder, Wn, 'bandpass');
            catch ME
                warning('Could not design Butterworth filter: %s. Skipping filtering.', ME.message);
                filteredSignals = processedSignals;
                return;
            end

            for ch = 1:nChannels
                x_ch = squeeze(processedSignals(:, :, ch)); % Samples x Windows
                % Filter along the time dimension (samples) for each window independently
                % Note: filtfilt processes columns by default. Need to process each window.
                % This assumes stationarity within the 1s window for filtering.
                % Alternative: filter the concatenated signal (might introduce edge artifacts between windows)
                
                % Filter each window (column)
                filtered_x_ch = zeros(size(x_ch));
                for win = 1:size(x_ch, 2)
                   filtered_x_ch(:, win) = filtfilt(b, a, x_ch(:, win));
                end
                filteredSignals(:, :, ch) = filtered_x_ch;

                 % Alternative: Transpose, filter rows, transpose back (might be faster?)
                 % filteredSignals(:,:,ch) = filtfilt(b, a, x_ch')';
            end
            fprintf('[%s] Antunes Filtering complete.\n', datetime);
        end

        function processedSignals = processSignals(obj, rawSignals, maxWindows)
             % Convenience method to apply both Zanoteli and Antunes steps.
             arguments
                obj
                rawSignals (:,:,:) {mustBeNumeric}
                maxWindows (1,1) {mustBeInteger, mustBePositive} = size(rawSignals, 2)
            end
            
            % 1. Apply Zanoteli
            zanoteliProcessed = obj.applyZanoteliPreprocessing(rawSignals, maxWindows);
            
            % 2. Apply Filtering
            processedSignals = obj.applyAntunesFiltering(zanoteliProcessed);
        end

    end
end