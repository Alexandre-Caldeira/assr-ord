classdef MetricsCalculator < handle
    % Calculates and visualizes performance metrics from pre-aggregated test results.

    properties (SetAccess = protected)
        paramList % Array of structs detailing the parameters for each aggregated row
        resultsTable % Table containing parameters and calculated rates/times per row
        nSignalFreqs % Number of signal frequencies used in tests
        nNoiseFreqs  % Number of noise frequencies used in tests
        nChannels    % Number of channels used (can be avg/typical if varies)

        avg_tp_rate, avg_fn_rate, avg_fp_rate, avg_tn_rate
        avg_min_detect_time, avg_avg_detect_time
        confmat_avg % Table summarizing average rates/times
    end

    methods
        function obj = MetricsCalculator(paramList, aggregatedData, nSignalFreqs, nNoiseFreqs, nChannels)
            % Constructor
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
            if isempty(paramList) || ~isstruct(aggregatedData), error('MetricsCalculator:InvalidInput', 'Input paramList or aggregatedData is empty or invalid.'); end
            obj.calculateRatesAndTimes(paramList, aggregatedData);
            obj.calculateAverageMetrics();
        end

        function calculateRatesAndTimes(obj, paramList, aggregatedData)
            % Calculates rates and average times per aggregated parameter set row.
            fprintf('[%s] Calculating final rates and average times...\n', datetime);
            nParamSets = numel(paramList);
            requiredFields = {'sumTP', 'sumFN', 'sumFP', 'sumTN', 'totalSignalOpps', 'totalNoiseOpps'};
            if ~all(isfield(aggregatedData, requiredFields)), error('MetricsCalculator:InvalidInput', 'aggregatedData missing required fields.'); end
            if numel(aggregatedData.sumTP) ~= nParamSets, error('MetricsCalculator:InvalidInput', 'Size mismatch between paramList and aggregatedData.'); end

            tp_rate = aggregatedData.sumTP ./ aggregatedData.totalSignalOpps;
            fn_rate = aggregatedData.sumFN ./ aggregatedData.totalSignalOpps;
            fp_rate = aggregatedData.sumFP ./ aggregatedData.totalNoiseOpps;
            tn_rate = aggregatedData.sumTN ./ aggregatedData.totalNoiseOpps;
            tp_rate(aggregatedData.totalSignalOpps == 0) = NaN; fn_rate(aggregatedData.totalSignalOpps == 0) = NaN;
            fp_rate(aggregatedData.totalNoiseOpps == 0) = NaN; tn_rate(aggregatedData.totalNoiseOpps == 0) = NaN;

            obj.resultsTable = struct2table(paramList); % Includes M, sim/exp params now
            obj.resultsTable.TP_Rate = tp_rate * 100; obj.resultsTable.FN_Rate = fn_rate * 100;
            obj.resultsTable.FP_Rate = fp_rate * 100; obj.resultsTable.TN_Rate = tn_rate * 100;

            if isfield(aggregatedData, 'sumMinDetTime') && isfield(aggregatedData, 'nMinDetTime')
                 avgMinTime = aggregatedData.sumMinDetTime ./ aggregatedData.nMinDetTime;
                 avgMinTime(aggregatedData.nMinDetTime == 0) = NaN;
                 obj.resultsTable.AvgMinDetectTime_s = avgMinTime;
            end
             if isfield(aggregatedData, 'sumAvgDetTime') && isfield(aggregatedData, 'nAvgDetTime')
                 avgAvgTime = aggregatedData.sumAvgDetTime ./ aggregatedData.nAvgDetTime;
                 avgAvgTime(aggregatedData.nAvgDetTime == 0) = NaN;
                 obj.resultsTable.AvgDetectTime_s = avgAvgTime;
            end
            fprintf('[%s] Rates and average times calculation complete.\n', datetime);
        end

        function calculateAverageMetrics(obj)
            % Calculates average metrics across all parameter sets in resultsTable.
             if isempty(obj.resultsTable), warning('Results table is empty.'); return; end
             % Average the rates calculated per parameter set (which might represent different groups)
             obj.avg_tp_rate = mean(obj.resultsTable.TP_Rate, 'omitnan');
             obj.avg_fn_rate = mean(obj.resultsTable.FN_Rate, 'omitnan');
             obj.avg_fp_rate = mean(obj.resultsTable.FP_Rate, 'omitnan');
             obj.avg_tn_rate = mean(obj.resultsTable.TN_Rate, 'omitnan');
             avg_vars = {'Avg FN Rate (%)', 'Avg FP Rate (%)', 'Avg TP Rate (%)', 'Avg TN Rate (%)'};
             avg_vals = [obj.avg_fn_rate, obj.avg_fp_rate, obj.avg_tp_rate, obj.avg_tn_rate];
             if ismember('AvgMinDetectTime_s', obj.resultsTable.Properties.VariableNames)
                 obj.avg_min_detect_time = mean(obj.resultsTable.AvgMinDetectTime_s, 'omitnan');
                 avg_vars = [avg_vars, 'Overall Avg Min Detect Time (s)']; avg_vals = [avg_vals, obj.avg_min_detect_time];
             else, obj.avg_min_detect_time = NaN; end
              if ismember('AvgDetectTime_s', obj.resultsTable.Properties.VariableNames)
                 obj.avg_avg_detect_time = mean(obj.resultsTable.AvgDetectTime_s, 'omitnan');
                 avg_vars = [avg_vars, 'Overall Avg Exam Time (TP only, s)']; avg_vals = [avg_vals, obj.avg_avg_detect_time];
              else, obj.avg_avg_detect_time = NaN; end
             obj.confmat_avg = array2table(avg_vals, 'VariableNames', avg_vars);
        end

        function displaySummary(obj, desiredAlpha)
            % Displays summary and validates FP rate.
             arguments, obj, desiredAlpha = [], end
             if isempty(obj.resultsTable), warning('No results to display.'); return; end
             disp('--- Results Summary (First 10 Rows) ---'); disp(head(obj.resultsTable, 10));
             disp('--- Overall Average Metrics Across Parameter Sets ---'); disp(obj.confmat_avg);
             if ~isempty(desiredAlpha) && ismember('FP_Rate', obj.resultsTable.Properties.VariableNames) && ~isnan(obj.avg_fp_rate)
                 fprintf('Validation: Desired Alpha = %.4f, Calculated Avg FP Rate = %.4f (%.2f%%)\n', desiredAlpha, obj.avg_fp_rate / 100, obj.avg_fp_rate);
                 % Check FP rate per group if grouping variable exists
                 groupingVars = {'noiseMean', 'subject', 'stimulus'}; % Potential grouping vars
                 groupVarFound = '';
                 for gv = groupingVars
                     if ismember(gv{1}, obj.resultsTable.Properties.VariableNames)
                         groupVarFound = gv{1};
                         break;
                     end
                 end
                 if ~isempty(groupVarFound)
                     try
                        groupStats = groupsummary(obj.resultsTable, groupVarFound, 'mean', 'FP_Rate');
                        fprintf('Average FP Rate per %s group:\n', groupVarFound);
                        disp(groupStats(:, {groupVarFound, 'mean_FP_Rate'}));
                        % Check if *each group's* avg FP rate is close
                        fp_dev = abs(groupStats.mean_FP_Rate / 100 - desiredAlpha);
                        if any(fp_dev > 0.015) % Increased tolerance
                           warning('FP rate in one or more groups deviates significantly from desired alpha!');
                        end
                     catch ME_group
                         warning('Could not calculate grouped FP rate summary: %s', ME_group.message);
                     end
                 elseif abs(obj.avg_fp_rate / 100 - desiredAlpha) > 0.015 % Check overall if no grouping found
                     warning('Calculated overall average FP rate deviates significantly from desired alpha!');
                 elseif ~isnan(obj.avg_fp_rate)
                     fprintf(' -> Overall FP rate matches desired alpha within tolerance.\n');
                 end

             elseif ~isempty(desiredAlpha) && ismember('FP_Rate', obj.resultsTable.Properties.VariableNames) && isnan(obj.avg_fp_rate) && obj.nNoiseFreqs > 0, warning('FP Rate is NaN, cannot validate against desired alpha (%.4f).', desiredAlpha); end
             fprintf('-------------------------------------------\n');
        end

        % --- Boxplot Methods ---
        function plotRateBoxcharts(obj, figHandle)
             if isempty(obj.resultsTable), warning('No results to plot.'); return; end
             if nargin < 2 || isempty(figHandle), figHandle = figure; else, figure(figHandle); end; clf(figHandle);
             subplot(2, 2, 1); boxchart(obj.resultsTable.TP_Rate); title('TP Rate (%)'); grid on; ylim([0 100]); ylabel('Rate (%)');
             subplot(2, 2, 2); boxchart(obj.resultsTable.FP_Rate); title('FP Rate (%)'); grid on; ylim([0 100]);
             subplot(2, 2, 3); boxchart(obj.resultsTable.FN_Rate); title('FN Rate (%)'); grid on; ylim([0 100]); ylabel('Rate (%)');
             subplot(2, 2, 4); boxchart(obj.resultsTable.TN_Rate); title('TN Rate (%)'); grid on; ylim([0 100]);
             sgtitle('Distribution of Rates Across All Parameter Sets (Boxcharts)');
        end

        function plotTimeBoxcharts(obj, figHandle)
             if isempty(obj.resultsTable), warning('No time results to plot.'); return; end
             if nargin < 2 || isempty(figHandle), figHandle = figure; else, figure(figHandle); end; clf(figHandle);
             hasMin = ismember('AvgMinDetectTime_s', obj.resultsTable.Properties.VariableNames);
             hasAvg = ismember('AvgDetectTime_s', obj.resultsTable.Properties.VariableNames);
             if ~hasMin && ~hasAvg, warning('No time results found in table.'); return; end
             if hasMin && hasAvg
                 subplot(1, 2, 1); boxchart(obj.resultsTable.AvgMinDetectTime_s); title('Avg Min Detection Time (s)'); grid on; ylabel('Time (s)');
                 subplot(1, 2, 2); boxchart(obj.resultsTable.AvgDetectTime_s); title('Avg Exam Time (TP only, s)'); grid on;
             elseif hasMin, boxchart(obj.resultsTable.AvgMinDetectTime_s); title('Avg Min Detection Time (s)'); grid on; ylabel('Time (s)');
             elseif hasAvg, boxchart(obj.resultsTable.AvgDetectTime_s); title('Avg Exam Time (TP only, s)'); grid on; ylabel('Time (s)');
             end
             sgtitle('Distribution of Average Detection Times Across Parameter Sets (Boxcharts)');
        end

        % --- Pareto Plotting Methods ---
        function plotParetoFront(obj, xMetric, yMetric, figHandle, options)
            % Plots Pareto front with optional labels and M-coloring.
             arguments
                 obj, xMetric char, yMetric char, figHandle = figure
                 options.addLabels logical = true
                 options.colorMetric char = 'M'
                 options.colormapName char = 'jet'
                 options.markerSizeData (1,1) {mustBePositive} = 20 % Size for background data points
                 options.markerSizeFront (1,1) {mustBePositive} = 50 % Size for pareto front points
                 options.markerAlphaData (1,1) {mustBeInRange(options.markerAlphaData, 0, 1)} = 0.7 % Alpha for background data points
             end
             if isempty(obj.resultsTable), warning('Cannot plot Pareto front. Table empty.'); return; end
             if ~ismember(xMetric, obj.resultsTable.Properties.VariableNames) || ~ismember(yMetric, obj.resultsTable.Properties.VariableNames)
                 warning('Cannot plot Pareto front. Metrics "%s", "%s" not found.', xMetric, yMetric); return;
             end
             useColor = ~isempty(options.colorMetric) && ismember(options.colorMetric, obj.resultsTable.Properties.VariableNames);
             if ~isempty(options.colorMetric) && ~useColor, warning('Color metric "%s" not found. Plotting without color.', options.colorMetric); end

             if nargin < 4 || isempty(figHandle), figHandle = figure; else, figure(figHandle); end; clf(figHandle); ax = gca;

             [paretoX, paretoY, xDataValid, yDataValid, originalFrontIndices] = obj.calculateParetoData(xMetric, yMetric);

             if isempty(paretoX)
                 if ~isempty(xDataValid)
                    scatter(ax, xDataValid, yDataValid, options.markerSizeData, 'k', 'filled', 'MarkerFaceAlpha', options.markerAlphaData, 'DisplayName', 'All Valid Parameter Sets');
                 end
                 xlabel(strrep(xMetric,'_',' ')); ylabel(strrep(yMetric,'_',' ')); title(sprintf('Pareto Front Calculation Failed/Empty: %s vs. %s', strrep(yMetric,'_',' '), strrep(xMetric,'_',' ')));
                 legend(ax, 'Location', 'best'); grid(ax, 'on'); return
             end

             colorDataValid = []; colorDataPareto = [];
             if useColor
                 colorVals = obj.resultsTable.(options.colorMetric);
                 colorDataValid = colorVals(obj.getValidIndices(xMetric, yMetric));
                 colorDataPareto = colorVals(originalFrontIndices);
             end

             hold(ax, 'on');
             if useColor && ~isempty(colorDataValid)
                 scatter(ax, xDataValid, yDataValid, options.markerSizeData, colorDataValid, 'filled', 'MarkerFaceAlpha', options.markerAlphaData, 'DisplayName', 'All Valid Parameter Sets');
             else
                 plot(ax, xDataValid, yDataValid, 'k.', 'MarkerSize', options.markerSizeData * 0.5, 'DisplayName', 'All Valid Parameter Sets'); % Smaller marker if no color
             end

             [sortedX, sort_idx] = sort(paretoX); sortedY = paretoY(sort_idx);
             sortedOriginalIndices = originalFrontIndices(sort_idx);
             if useColor && ~isempty(colorDataPareto), sortedColorDataPareto = colorDataPareto(sort_idx); end

             plot(ax, sortedX, sortedY, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'HandleVisibility', 'off');
             if useColor && ~isempty(sortedColorDataPareto)
                 scatter(ax, sortedX, sortedY, options.markerSizeFront, sortedColorDataPareto, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'Pareto Front');
             else
                 plot(ax, sortedX, sortedY, 'ob', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0, 'MarkerSize', sqrt(options.markerSizeFront), 'DisplayName', 'Pareto Front'); % Adjust marker size
             end

             if useColor && ~isempty(colorDataValid)
                 try, colormap(ax, options.colormapName); catch, colormap(ax, 'jet'); warning('Colormap "%s" invalid, using jet.', options.colormapName); end % Use try-catch for colormap
                 cb = colorbar(ax); cb.Label.String = strrep(options.colorMetric,'_',' ');
                 if ~isempty(colorDataValid), caxis(ax, [min(colorDataValid), max(colorDataValid)]); end
             end

             if options.addLabels
                 frontTable = obj.resultsTable(sortedOriginalIndices,:);
                 for i = 1:numel(sortedX)
                     row = frontTable(i,:); labelParts = {};
                     if ismember('K', row.Properties.VariableNames), labelParts{end+1} = sprintf('K=%d', row.K); end
                     if ismember('startWindow', row.Properties.VariableNames), labelParts{end+1} = sprintf('SW=%d', row.startWindow); end
                     if ismember('M', row.Properties.VariableNames), labelParts{end+1} = sprintf('M=%d', row.M); end
                     if ismember('noiseMean', row.Properties.VariableNames), labelParts{end+1} = sprintf('nM=%.0f', row.noiseMean); end
                     % if ismember('noiseStd', row.Properties.VariableNames), labelParts{end+1} = sprintf('nS=%.0f', row.noiseStd); end
                     % if ismember('subject', row.Properties.VariableNames), labelParts{end+1} = sprintf('S=%d', row.subject); end
                     % if ismember('stimulus', row.Properties.VariableNames), labelParts{end+1} = sprintf('St=%d', row.stimulus); end
                     labelString = strjoin(labelParts, ',');
                     text(ax, sortedX(i), sortedY(i), [' ' labelString], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 7, 'Interpreter', 'none');
                 end
             end

             hold(ax, 'off'); xlabel(strrep(xMetric,'_',' ')); ylabel(strrep(yMetric,'_',' '));
             title(sprintf('Pareto Front: %s vs. %s', strrep(yMetric,'_',' '), strrep(xMetric,'_',' ')));
             legend(ax, 'Location', 'best', 'Interpreter', 'none'); grid(ax, 'on');
        end

        function plotGroupedParetoFront(obj, groupingVar, xMetric, yMetric, figHandle, options)
            % Plots multiple Pareto fronts, colored by M, grouped by groupingVar.
             arguments
                 obj, groupingVar char, xMetric char, yMetric char, figHandle = figure
                 options.colormapName char = 'jet'
                 options.markerSizeFront (1,1) {mustBePositive} = 36 % Size for pareto front points
             end
             if isempty(obj.resultsTable), warning('Cannot plot grouped Pareto. Table empty.'); return; end
             if ~ismember(groupingVar, obj.resultsTable.Properties.VariableNames), warning('Grouping variable "%s" not found.', groupingVar); return; end
             if ~ismember(xMetric, obj.resultsTable.Properties.VariableNames) || ~ismember(yMetric, obj.resultsTable.Properties.VariableNames)
                 warning('Cannot plot grouped Pareto. Metrics "%s", "%s" not found.', xMetric, yMetric); return; end
             if ~ismember('M', obj.resultsTable.Properties.VariableNames), warning('Cannot color by M, column not found.'); useColor = false; else, useColor = true; end

             if nargin < 5 || isempty(figHandle), figHandle = figure; else, figure(figHandle); end; clf(figHandle); ax = gca;

             groupValues = unique(obj.resultsTable.(groupingVar)); numGroups = numel(groupValues);
             groupLineColors = lines(numGroups); % Colors for lines distinguishing groups
             hold(ax, 'on'); legendEntries = cell(numGroups, 1); allColorData = [];

             for i = 1:numGroups
                 groupVal = groupValues(i); groupMask = obj.resultsTable.(groupingVar) == groupVal;
                 groupTable = obj.resultsTable(groupMask,:);
                 [paretoX, paretoY, ~, ~, originalFrontIndices] = obj.calculateParetoData(xMetric, yMetric, groupTable);

                 if ~isempty(paretoX)
                     [sortedX, sort_idx] = sort(paretoX); sortedY = paretoY(sort_idx);
                     sortedOriginalIndicesGroup = originalFrontIndices(sort_idx);
                     plot(ax, sortedX, sortedY, '-', 'LineWidth', 1.0, 'Color', groupLineColors(i,:), 'HandleVisibility', 'off');

                     colorDataParetoGroup = [];
                     if useColor
                         colorDataParetoGroup = groupTable.M(sortedOriginalIndicesGroup);
                         allColorData = [allColorData; colorDataParetoGroup];
                         scatter(ax, sortedX, sortedY, options.markerSizeFront, colorDataParetoGroup, 'filled', ...
                                 'MarkerEdgeColor', groupLineColors(i,:)*0.7, 'LineWidth', 0.7, ... % Darker edge based on group color
                                 'DisplayName', sprintf('%s = %g', groupingVar, groupVal));
                     else
                         plot(ax, sortedX, sortedY, 'o', 'MarkerSize', sqrt(options.markerSizeFront), ...
                              'MarkerFaceColor', groupLineColors(i,:), 'MarkerEdgeColor', 'k', ...
                              'DisplayName', sprintf('%s = %g', groupingVar, groupVal));
                     end
                     legendEntries{i} = sprintf('%s = %g', groupingVar, groupVal);
                 else
                     warning('No Pareto front found for group %s = %g.', groupingVar, groupVal);
                     legendEntries{i} = sprintf('%s = %g (No Front)', groupingVar, groupVal);
                 end
             end

             if useColor && ~isempty(allColorData)
                 try, colormap(ax, options.colormapName); catch, colormap(ax, 'jet'); warning('Colormap "%s" invalid, using jet.', options.colormapName); end
                 cb = colorbar(ax); cb.Label.String = 'Window Size (M)';
                 if ~isempty(allColorData), caxis(ax, [min(allColorData), max(allColorData)]); end
             end

             hold(ax, 'off'); xlabel(strrep(xMetric,'_',' ')); ylabel(strrep(yMetric,'_',' '));
             title(sprintf('Pareto Fronts: %s vs. %s (Grouped by %s)', strrep(yMetric,'_',' '), strrep(xMetric,'_',' '), groupingVar));
             legend(legendEntries(~cellfun('isempty',legendEntries)), 'Location', 'best', 'Interpreter', 'none'); grid(ax, 'on');
        end

        function plotAllMetrics(obj, options)
            % Generates a standard set of plots using Name-Value pairs.
            arguments
                obj
                options.includeBoxcharts logical = true
                options.includeParetoFPTP logical = true
                options.includeParetoTimeTP logical = true
                options.includeGroupedPareto logical = false
                options.groupingVar char = 'noiseMean'
                options.paretoColormap char = 'jet'
                options.addParetoLabels logical = true
            end

            fprintf('Generating standard plots...\n');
            existingFigs = findall(groot, 'Type', 'figure');
            if isempty(existingFigs), baseFigNum = 1;
            else, figNumbers = [existingFigs.Number]; figNumbers = figNumbers(floor(figNumbers) == figNumbers & ~isempty(figNumbers));
                 if isempty(figNumbers), baseFigNum = 1; else, baseFigNum = max(figNumbers) + 1; end
            end
            figCount = 0;

            if options.includeBoxcharts
                obj.plotRateBoxcharts(figure(baseFigNum + figCount)); figCount = figCount + 1;
                if ismember('AvgMinDetectTime_s', obj.resultsTable.Properties.VariableNames)
                    obj.plotTimeBoxcharts(figure(baseFigNum + figCount)); figCount = figCount + 1;
                end
            end
             if options.includeParetoFPTP
                 obj.plotParetoFront('FP_Rate', 'TP_Rate', figure(baseFigNum + figCount), ...
                     'addLabels', options.addParetoLabels, 'colormapName', options.paretoColormap);
                 figCount = figCount + 1;
             end
             if options.includeParetoTimeTP && ismember('AvgDetectTime_s', obj.resultsTable.Properties.VariableNames)
                 obj.plotParetoFront('TP_Rate', 'AvgDetectTime_s', figure(baseFigNum + figCount), ... % TP (X) vs Time (Y)
                      'addLabels', options.addParetoLabels, 'colormapName', options.paretoColormap);
                 figCount = figCount + 1;
             end
             if options.includeGroupedPareto
                 if ismember(options.groupingVar, obj.resultsTable.Properties.VariableNames) && ismember('AvgDetectTime_s', obj.resultsTable.Properties.VariableNames)
                     obj.plotGroupedParetoFront(options.groupingVar, 'TP_Rate', 'AvgDetectTime_s', figure(baseFigNum + figCount), ... % TP(X) vs Time(Y) grouped
                          'colormapName', options.paretoColormap);
                     figCount = figCount + 1;
                 else
                     warning('Cannot generate grouped Pareto plot: Grouping variable "%s" or metric "AvgDetectTime_s" not found.', options.groupingVar);
                 end
             end
             fprintf('Plot generation finished.\n');
        end

        function results = getResultsTable(obj), results = obj.resultsTable; end
        function avgResults = getAverageResults(obj), avgResults = obj.confmat_avg; end
    end

    methods (Access = private)
        function [paretoX, paretoY, xDataValid, yDataValid, originalFrontIndices] = calculateParetoData(obj, xMetric, yMetric, dataTable)
            % Calculates Pareto front data points and indices.
            arguments, obj, xMetric char, yMetric char, dataTable = obj.resultsTable, end
            paretoX = []; paretoY = []; xDataValid = []; yDataValid = []; originalFrontIndices = [];
            if isempty(dataTable), return; end

            minimizeX = true; if contains(xMetric,'TP_Rate','IgnoreCase',true), minimizeX=false; end
            maximizeY = true; if contains(yMetric,'FP_Rate','IgnoreCase',true) || contains(yMetric,'FN_Rate','IgnoreCase',true) || contains(yMetric,'Time','IgnoreCase',true), maximizeY=false; end

            xVals = dataTable.(xMetric); yVals = dataTable.(yMetric);
            valid_idx_logical = ~isnan(xVals) & ~isnan(yVals);
            valid_idx_numeric = find(valid_idx_logical);

            if sum(valid_idx_logical) < 2, return; end

            xDataValid = xVals(valid_idx_logical); yDataValid = yVals(valid_idx_logical);
            obj1 = xDataValid; if ~minimizeX, obj1 = -obj1; end
            obj2 = yDataValid; if maximizeY, obj2 = -obj2; end

            try
                [frontIndicesRel, ~] = paretoFront([obj1, obj2]);
                if isempty(frontIndicesRel), return; end
                originalFrontIndices = valid_idx_numeric(frontIndicesRel); % Indices relative to dataTable
                paretoX = xDataValid(frontIndicesRel); paretoY = yDataValid(frontIndicesRel);
            catch ME_pareto
                warning('MetricsCalculator:ParetoError', 'Pareto calculation failed for %s vs %s: %s', xMetric, yMetric, ME_pareto.message);
            end
        end

        function idx = getValidIndices(obj, xMetric, yMetric, dataTable)
             arguments, obj, xMetric char, yMetric char, dataTable = obj.resultsTable, end
             if isempty(dataTable) || ~ismember(xMetric, dataTable.Properties.VariableNames) || ~ismember(yMetric, dataTable.Properties.VariableNames)
                 idx = false(height(dataTable), 1); return;
             end
             xVals = dataTable.(xMetric); yVals = dataTable.(yMetric);
             idx = ~isnan(xVals) & ~isnan(yVals);
        end
    end
end