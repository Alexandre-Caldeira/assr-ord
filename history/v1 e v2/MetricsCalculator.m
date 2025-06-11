classdef MetricsCalculator < handle
    % Calculates and visualizes performance metrics from aggregated test results.

    properties (SetAccess = protected)
        resultsTable % Table containing parameters and aggregated counts/rates per configuration
        nSignalFreqs % Number of signal frequencies used in tests
        nNoiseFreqs  % Number of noise frequencies used in tests
        nChannels    % Number of channels used in tests

        % Calculated overall metrics
        avg_tp_rate
        avg_fn_rate
        avg_fp_rate
        avg_tn_rate
        avg_min_detect_time
        avg_avg_detect_time
        confmat_avg % Table summarizing average rates/times
    end

    methods
        function obj = MetricsCalculator(paramList, aggregatedData, nSignalFreqs, nNoiseFreqs, nChannels)
            % Constructor for MetricsCalculator
            % paramList: Array of structs containing parameters for each configuration
            % aggregatedData: Struct containing aggregated counts/times vectors
            %   - TP, FN, FP, TN (counts per config)
            %   - (Optional) MinDetectTime, AvgDetectTime (per config)
            %   - (Optional) nSignalOpps, nNoiseOpps (opportunities per config)
            % nSignalFreqs, nNoiseFreqs, nChannels: Constants used during testing

            arguments
                paramList struct
                aggregatedData struct
                nSignalFreqs (1,1) {mustBeInteger, mustBeNonnegative}
                nNoiseFreqs (1,1) {mustBeInteger, mustBeNonnegative}
                nChannels (1,1) {mustBeInteger, mustBePositive}
            end

            obj.nSignalFreqs = nSignalFreqs;
            obj.nNoiseFreqs = nNoiseFreqs;
            obj.nChannels = nChannels;

            if isempty(paramList) || ~isstruct(aggregatedData)
                error('MetricsCalculator:InvalidInput', 'Input paramList or aggregatedData is empty or invalid.');
            end

            obj.calculateMetrics(paramList, aggregatedData);
            obj.calculateAverageMetrics();
        end

        function calculateMetrics(obj, paramList, aggregatedData)
            % Calculates rates per configuration and populates resultsTable.
            fprintf('[%s] Calculating performance metrics...\n', datetime);

            nConfigs = numel(paramList);
            requiredFields = {'TP', 'FN', 'FP', 'TN'};
            if ~all(isfield(aggregatedData, requiredFields))
                 error('MetricsCalculator:InvalidInput', 'aggregatedData missing required fields (TP, FN, FP, TN).');
            end
            
            % --- Calculate Denominators (Opportunities) ---
            % Check if provided, otherwise calculate
            if isfield(aggregatedData, 'nSignalOpps')
                nSignalOpps = aggregatedData.nSignalOpps;
                if numel(nSignalOpps) ~= nConfigs
                    error('MetricsCalculator:InvalidInput', 'Size mismatch for nSignalOpps.');
                end
            else
                nSignalOpps = repmat(obj.nSignalFreqs * obj.nChannels, nConfigs, 1);
            end
            
            if isfield(aggregatedData, 'nNoiseOpps')
                nNoiseOpps = aggregatedData.nNoiseOpps;
                 if numel(nNoiseOpps) ~= nConfigs
                    error('MetricsCalculator:InvalidInput', 'Size mismatch for nNoiseOpps.');
                 end
            else
                nNoiseOpps = repmat(obj.nNoiseFreqs * obj.nChannels, nConfigs, 1);
            end


            % --- Calculate Rates ---
            tp_rate = aggregatedData.TP ./ nSignalOpps;
            fn_rate = aggregatedData.FN ./ nSignalOpps;
            fp_rate = aggregatedData.FP ./ nNoiseOpps;
            tn_rate = aggregatedData.TN ./ nNoiseOpps;

            % Handle division by zero
            tp_rate(nSignalOpps == 0) = NaN;
            fn_rate(nSignalOpps == 0) = NaN;
            fp_rate(nNoiseOpps == 0) = NaN;
            tn_rate(nNoiseOpps == 0) = NaN;

            % --- Create Table ---
            obj.resultsTable = struct2table(paramList);
            obj.resultsTable.TP_Rate = tp_rate * 100; % In percent
            obj.resultsTable.FN_Rate = fn_rate * 100;
            obj.resultsTable.FP_Rate = fp_rate * 100;
            obj.resultsTable.TN_Rate = tn_rate * 100;

            % Add optional time columns if present
            if isfield(aggregatedData, 'MinDetectTime')
                obj.resultsTable.MinDetectTime_s = aggregatedData.MinDetectTime;
            end
            if isfield(aggregatedData, 'AvgDetectTime')
                obj.resultsTable.AvgDetectTime_s = aggregatedData.AvgDetectTime;
            end

            fprintf('[%s] Metrics calculation complete. Results table created.\n', datetime);
        end

        function calculateAverageMetrics(obj)
            % Calculates average metrics across all configurations.
            if isempty(obj.resultsTable)
                warning('Results table is empty. Cannot calculate average metrics.');
                return;
            end

            obj.avg_tp_rate = mean(obj.resultsTable.TP_Rate, 'omitnan');
            obj.avg_fn_rate = mean(obj.resultsTable.FN_Rate, 'omitnan');
            obj.avg_fp_rate = mean(obj.resultsTable.FP_Rate, 'omitnan');
            obj.avg_tn_rate = mean(obj.resultsTable.TN_Rate, 'omitnan');

            avg_vars = {'Avg FN Rate (%)', 'Avg FP Rate (%)', 'Avg TP Rate (%)', 'Avg TN Rate (%)'};
            avg_vals = [obj.avg_fn_rate, obj.avg_fp_rate, obj.avg_tp_rate, obj.avg_tn_rate];

            if ismember('MinDetectTime_s', obj.resultsTable.Properties.VariableNames)
                obj.avg_min_detect_time = mean(obj.resultsTable.MinDetectTime_s, 'omitnan');
                avg_vars = [avg_vars, 'Avg Min Detect Time (s)'];
                avg_vals = [avg_vals, obj.avg_min_detect_time];
            else
                 obj.avg_min_detect_time = NaN;
            end

             if ismember('AvgDetectTime_s', obj.resultsTable.Properties.VariableNames)
                obj.avg_avg_detect_time = mean(obj.resultsTable.AvgDetectTime_s, 'omitnan');
                avg_vars = [avg_vars, 'Avg Avg Detect Time (s)'];
                avg_vals = [avg_vals, obj.avg_avg_detect_time];
             else
                 obj.avg_avg_detect_time = NaN;
             end

            obj.confmat_avg = array2table(avg_vals, 'VariableNames', avg_vars);
        end

        function displaySummary(obj)
            % Displays the summary results table and average confusion matrix.
             if isempty(obj.resultsTable)
                warning('No results to display.');
                return;
            end
            disp('Results Summary Table (Rates in %):');
            disp(head(obj.resultsTable, 10)); % Show first 10 rows

            disp('Overall Average Confusion Rates (%) and Detection Times (s):');
            disp(obj.confmat_avg);
        end

        function plotRateDistributions(obj, figHandle)
            % Plots box charts of the rate distributions.
             if isempty(obj.resultsTable)
                warning('No results to plot.');
                return;
            end
            if nargin < 2 || isempty(figHandle)
                figHandle = figure;
            else
                figure(figHandle); % Bring to front or create if handle invalid
            end
            clf(figHandle); % Clear figure

            subplot(2, 2, 1); boxchart(obj.resultsTable.TP_Rate); title('TP Rate (%)'); grid on; ylim([0 100]);
            subplot(2, 2, 2); boxchart(obj.resultsTable.FP_Rate); title('FP Rate (%)'); grid on; ylim([0 100]);
            subplot(2, 2, 3); boxchart(obj.resultsTable.FN_Rate); title('FN Rate (%)'); grid on; ylim([0 100]);
            subplot(2, 2, 4); boxchart(obj.resultsTable.TN_Rate); title('TN Rate (%)'); grid on; ylim([0 100]);
            sgtitle('Distribution of Rates Across All Parameter Sets');
        end

        function plotTimeDistributions(obj, figHandle)
            % Plots box charts of the detection time distributions.
             if isempty(obj.resultsTable) || ...
                ~ismember('MinDetectTime_s', obj.resultsTable.Properties.VariableNames)
                warning('No detection time results to plot.');
                return;
            end
             if nargin < 2 || isempty(figHandle)
                figHandle = figure;
            else
                figure(figHandle);
            end
            clf(figHandle);

            hasMin = ismember('MinDetectTime_s', obj.resultsTable.Properties.VariableNames);
            hasAvg = ismember('AvgDetectTime_s', obj.resultsTable.Properties.VariableNames);

            if hasMin && hasAvg
                subplot(1, 2, 1); boxchart(obj.resultsTable.MinDetectTime_s); title('Min Detection Time (s)'); grid on;
                subplot(1, 2, 2); boxchart(obj.resultsTable.AvgDetectTime_s); title('Avg Detection Time (s)'); grid on;
            elseif hasMin
                 boxchart(obj.resultsTable.MinDetectTime_s); title('Min Detection Time (s)'); grid on;
            elseif hasAvg
                 boxchart(obj.resultsTable.AvgDetectTime_s); title('Avg Detection Time (s)'); grid on;
            end
            sgtitle('Distribution of Detection Times (TP only) Across All Parameter Sets');
        end

        function plotParetoFront(obj, xMetric, yMetric, figHandle)
            % Plots a Pareto front for two specified metrics from the results table.
            % Example: plotParetoFront('MinDetectTime_s', 'TP_Rate', figure)
             arguments
                 obj
                 xMetric char % Variable name in resultsTable for X-axis (e.g., 'MinDetectTime_s', 'FP_Rate')
                 yMetric char % Variable name in resultsTable for Y-axis (e.g., 'TP_Rate')
                 figHandle = figure % Optional figure handle
             end

             if isempty(obj.resultsTable) || ...
                ~ismember(xMetric, obj.resultsTable.Properties.VariableNames) || ...
                ~ismember(yMetric, obj.resultsTable.Properties.VariableNames)
                 warning('Cannot plot Pareto front. Table empty or metrics "%s", "%s" not found.', xMetric, yMetric);
                 return;
             end
             if nargin < 4 || isempty(figHandle)
                figHandle = figure;
             else
                figure(figHandle);
             end
             clf(figHandle);

             % Determine optimization direction (assume lower X and higher Y is better)
             % Adjust if needed (e.g., if minimizing TP rate)
             minimizeX = true;
             maximizeY = true; % Typical for TP Rate

             xVals = obj.resultsTable.(xMetric);
             yVals = obj.resultsTable.(yMetric);

             valid_idx = ~isnan(xVals) & ~isnan(yVals);
             if sum(valid_idx) < 2
                 warning('Not enough valid data points (need >= 2) to generate Pareto plot for %s vs %s.', xMetric, yMetric);
                 return;
             end

             xData = xVals(valid_idx);
             yData = yVals(valid_idx);

             % Prepare data for paretoFront (minimization)
             obj1 = xData;
             if ~minimizeX
                 obj1 = -obj1; % Minimize negative if we want to maximize X
             end

             obj2 = yData;
             if maximizeY
                 obj2 = -obj2; % Minimize negative if we want to maximize Y
             end

             try
                [p, p_indices_rel] = paretoFront([obj1, obj2]); % p contains minimized objectives
             catch ME_pareto
                 warning('Could not calculate Pareto front: %s', ME_pareto.message);
                 % Fallback: just plot the points
                 plot(xData, yData, 'k.', 'MarkerSize', 8, 'DisplayName', 'All Configurations');
                 xlabel(strrep(xMetric,'_',' ')); ylabel(strrep(yMetric,'_',' '));
                 title(sprintf('Scatter Plot: %s vs. %s', strrep(yMetric,'_',' '), strrep(xMetric,'_',' ')));
                 grid on; legend show;
                 return;
             end


             % Convert back to original scales for plotting
             paretoX = p(:, 1);
             if ~minimizeX
                 paretoX = -paretoX;
             end
             paretoY = p(:, 2);
             if maximizeY
                 paretoY = -paretoY;
             end

             % Sort pareto points for cleaner line plot (e.g., by X value)
             [sortedX, sort_idx] = sort(paretoX);
             sortedY = paretoY(sort_idx);

             % Plot all points and the Pareto front
             plot(xData, yData, 'k.', 'MarkerSize', 8, 'DisplayName', 'All Configurations');
             hold on;
             plot(sortedX, sortedY, '-ob', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'Pareto Front');
             hold off;
             xlabel(strrep(xMetric,'_',' '));
             ylabel(strrep(yMetric,'_',' '));
             title(sprintf('Pareto Front: %s vs. %s', strrep(yMetric,'_',' '), strrep(xMetric,'_',' ')));
             legend('Location', 'best');
             grid on;
             % Optional: Adjust axis limits if needed
             % xlim([min(xData) max(xData)]); ylim([min(yData) max(yData)]);
        end

        function results = getResultsTable(obj)
            % Returns the calculated results table.
            results = obj.resultsTable;
        end

         function avgResults = getAverageResults(obj)
            % Returns the calculated average results table.
            avgResults = obj.confmat_avg;
        end

    end % End methods
end % End classdef