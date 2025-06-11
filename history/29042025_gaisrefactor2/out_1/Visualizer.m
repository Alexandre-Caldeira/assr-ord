% --- FILE START: Visualizer.m ---
classdef Visualizer < handle
    % Handles plotting of results from the sequential testing workflow.

    properties
        DataManager DataManager % Access to metadata (frequencies, mode, etc.)
        MSCAnalyzer MSCAnalyzer % Access to MSC values and epoch definitions
        SequentialTester SequentialTester % Access to test decisions, stopping stages, performance metrics
    end

    methods
        function obj = Visualizer(dataManager, mscAnalyzer, sequentialTester)
            % Constructor: Stores handles to the data and analysis objects.
            arguments
                dataManager DataManager
                mscAnalyzer MSCAnalyzer
                sequentialTester SequentialTester
            end
            obj.DataManager = dataManager;
            obj.MSCAnalyzer = mscAnalyzer;
            obj.SequentialTester = sequentialTester;
            fprintf('Visualizer: Initialized.\n');
        end

        function figHandle = plotMSCandThresholds(obj, idx1, idx2, p1Idx, p2Idx, channelIdx)
            % Plots MSC accumulation vs. thresholds for one specific exam/param set/channel.
            % Args:
            %   idx1: Stimulus/Mean index.
            %   idx2: Subject/Std index.
            %   p1Idx: Index for the first parameter varied in MSCAnalyzer.
            %   p2Idx: Index for the second parameter varied in MSCAnalyzer.
            %   channelIdx: Index of the channel to plot.
            % Returns:
            %   figHandle: Handle to the created figure. Returns empty if plotting fails.

            figHandle = []; % Initialize output

            % --- Input Validation and Data Retrieval ---
            if isempty(obj.SequentialTester) || isempty(obj.SequentialTester.TestDecisions) || isempty(obj.SequentialTester.StrategyParamsArray)
                warning('Visualizer: SequentialTester results or parameters are missing. Run the test first.'); return;
            end
             if isempty(obj.MSCAnalyzer) || isempty(obj.MSCAnalyzer.mscValues)
                warning('Visualizer: MSCAnalyzer data is missing. Run MSC analysis first.'); return;
            end


            mscData = obj.MSCAnalyzer.getMSCData(idx1, idx2, p1Idx, p2Idx);
            decisions = obj.SequentialTester.TestDecisions{idx1, idx2, p1Idx, p2Idx};
            stopping = obj.SequentialTester.StoppingStages{idx1, idx2, p1Idx, p2Idx};
            strategyParams = obj.SequentialTester.StrategyParamsArray{idx1, idx2, p1Idx, p2Idx};
            epochs = obj.MSCAnalyzer.getEpochs(idx1, idx2, p1Idx, p2Idx);
            M_vals = obj.MSCAnalyzer.getWindowsPerStage(idx1, idx2, p1Idx, p2Idx);

            if isempty(mscData) || isempty(decisions) || isempty(stopping) || isempty(epochs) || isempty(M_vals) || isempty(strategyParams)
                 warning('Visualizer: Data missing for the specified indices (Idx1:%d, Idx2:%d, P1:%d, P2:%d). Cannot plot.', idx1, idx2, p1Idx, p2Idx); return;
            end

            [nBins, K, nChan] = size(mscData);
            if channelIdx < 1 || channelIdx > nChan
                error('Visualizer: Selected channel index %d is out of bounds [1, %d].', channelIdx, nChan);
            end
            if K == 0; warning('Visualizer: MSC data has 0 stages. Cannot plot.'); return; end

            % --- Recalculate Thresholds for this specific test ---
            analysisInfo.K = K;
            analysisInfo.M = M_vals;
            analysisInfo.MM = epochs(2:end) - 1; % Required for ETS
            try
                thresholds = obj.SequentialTester.Strategy.calculateThresholds(strategyParams, analysisInfo);
            catch ME
                 warning('Visualizer: Could not recalculate thresholds for plotting. Error: %s', ME.message); return;
            end

            % --- Prepare Data for Plotting ---
            mscSequenceChannel = squeeze(mscData(:, :, channelIdx));
            cumulativeMSC = cumsum(mscSequenceChannel, 2); % Accumulate MSC for the selected channel

            % Identify signal/noise bins (assuming simulation mode)
            sigFreqBins = []; noiseFreqBins = [];
            if strcmp(obj.DataManager.Mode, 'sim')
                 fs = obj.DataManager.Fs;
                 nfft = obj.DataManager.getNFFT();
                 freqResolution = fs / nfft;
                 sigFreqs_Hz = obj.DataManager.SimSignalFrequencies;
                 noiseFreqs_Hz = obj.DataManager.SimNoiseFrequencies;
                 sigFreqBins = round(sigFreqs_Hz / freqResolution) + 1;
                 noiseFreqBins = round(noiseFreqs_Hz / freqResolution) + 1;
                 maxBin = obj.DataManager.getNumBins();
                 sigFreqBins = unique(sigFreqBins(sigFreqBins >= 1 & sigFreqBins <= maxBin));
                 noiseFreqBins = unique(noiseFreqBins(noiseFreqBins >= 1 & noiseFreqBins <= maxBin));
                 noiseFreqBins = setdiff(noiseFreqBins, sigFreqBins);
            else
                warning('Visualizer: Not in simulation mode. Cannot distinguish signal/noise frequencies for plotting colors.');
                % Plot all bins with a default style
                sigFreqBins = 1:nBins; % Treat all as 'signal' for plotting purposes
            end

            % --- Create Plot ---
            figHandle = figure; hold on; grid on; box on;
             param1Name = obj.MSCAnalyzer.analyzedParam1Name;
             param1Val = obj.MSCAnalyzer.analyzedParam1Values(p1Idx);
             param2Name = obj.MSCAnalyzer.analyzedParam2Name;
             param2Val = obj.MSCAnalyzer.analyzedParam2Values(p2Idx);

            plotTitle = sprintf('Cumulative MSC vs. Thresholds (%s = %.1f, %s = %.1f)\nIdx1=%d, Idx2=%d, Chan=%d, Strategy=%s', ...
                               param1Name, param1Val, param2Name, param2Val, ...
                               idx1, idx2, channelIdx, obj.SequentialTester.Strategy.Name);
            title(plotTitle);
            xlabel('Stage (k)');
             if strcmp(obj.SequentialTester.Strategy.Name,'CGST_Beta')
                 ylabel('Cumulative MSC'); % CGST uses cumulative
             else
                 ylabel('Per-Stage MSC'); % ETS uses per-stage
                 cumulativeMSC = mscSequenceChannel; % Overwrite for plotting ETS
             end

            hPlots = []; % Handles for legend
            legendEntries = {};

            % Plot Signal Frequencies
            if ~isempty(sigFreqBins)
                 validSigBins = sigFreqBins(sigFreqBins <= nBins); % Ensure bins are within actual data
                 if ~isempty(validSigBins)
                    hPlots(end+1) = plot(1:K, cumulativeMSC(validSigBins, :)', '-b'); % Plot all signal freqs in blue
                    legendEntries{end+1} = 'Signal Freq Bins';
                 end
            end

            % Plot Noise Frequencies
             if ~isempty(noiseFreqBins)
                 validNoiseBins = noiseFreqBins(noiseFreqBins <= nBins); % Ensure bins are within actual data
                  if ~isempty(validNoiseBins)
                    hPlots(end+1) = plot(1:K, cumulativeMSC(validNoiseBins, :)', ':r'); % Plot all noise freqs in red dotted
                    legendEntries{end+1} = 'Noise Freq Bins';
                  end
             end

            % Plot Thresholds
            if isfield(thresholds, 'efficacy') % For CGST
                hPlots(end+1) = plot(1:K, thresholds.efficacy, '--g', 'LineWidth', 2);
                legendEntries{end+1} = 'Efficacy Threshold';
            end
             if isfield(thresholds, 'detection') % For ETS
                hPlots(end+1) = plot(1:K, thresholds.detection, '--g', 'LineWidth', 2);
                legendEntries{end+1} = 'Detection Threshold';
                if isfield(thresholds,'NDC') && thresholds.NDC > 1
                    legendEntries{end} = sprintf('Detection Threshold (NDC=%d)', thresholds.NDC);
                end
            end

            if obj.SequentialTester.Strategy.RequiresFutilityThreshold && isfield(thresholds, 'futility') && ~isempty(thresholds.futility)
                hPlots(end+1) = plot(1:K, thresholds.futility, '--m', 'LineWidth', 2);
                legendEntries{end+1} = 'Futility Threshold';
            end

            % Add legend (only use unique handles if lines were plotted multiple times)
            if ~isempty(hPlots)
                [~, uniqueIdx] = unique(hPlots); % Get indices of unique plot handles
                 legend(hPlots(uniqueIdx), legendEntries(uniqueIdx));
            end

            % Adjust Y-axis limits
            allDataToPlot = cumulativeMSC(:);
             if isfield(thresholds, 'efficacy'); allDataToPlot = [allDataToPlot; thresholds.efficacy(:)]; end
             if isfield(thresholds, 'detection'); allDataToPlot = [allDataToPlot; thresholds.detection(:)]; end
             if isfield(thresholds, 'futility') && ~isempty(thresholds.futility); allDataToPlot = [allDataToPlot; thresholds.futility(:)]; end
             validData = allDataToPlot(~isnan(allDataToPlot) & ~isinf(allDataToPlot));
             if ~isempty(validData)
                maxY = max(validData);
                minY = min(validData);
                ylim([min(0, minY - 0.1*abs(minY)), max(maxY + 0.1*abs(maxY), 0.1)]); % Dynamic ylim
             else
                 ylim([0, 1]); % Default if no valid data
             end
             xlim([1, K]);

            hold off;
        end

        function figHandle = plotFalsePositiveRate(obj, p)
            % Plots False Positive Rate (FPR) across the parameter sets tested.
            % Requires performance metrics to be calculated first (simulation mode).
            arguments
                obj
                p.PlotType {mustBeMember(p.PlotType, {'line', 'image'})} = 'image'
            end

            figHandle = [];

             if isempty(obj.SequentialTester) || isempty(obj.SequentialTester.PerformanceMetrics) || ~isfield(obj.SequentialTester.PerformanceMetrics, 'PerParameter')
                 warning('Visualizer: Performance metrics not calculated. Cannot plot FPR. Run calculatePerformance first.'); return;
             end
             paramResults = obj.SequentialTester.PerformanceMetrics.PerParameter;
             [nParam1, nParam2] = size(paramResults);
             if nParam1 == 0 || nParam2 == 0; warning('Visualizer: No parameter results found.'); return; end

             fpr_matrix = zeros(nParam1, nParam2) * NaN; % Initialize with NaN
             for i = 1:nParam1
                 for j = 1:nParam2
                     % Check if TotalNoiseOpportunities is present and positive
                     if isfield(paramResults(i,j), 'TotalNoiseOpportunities') && paramResults(i,j).TotalNoiseOpportunities > 0
                         % Check if FP field exists
                         if isfield(paramResults(i,j), 'FP')
                             fpr_matrix(i,j) = 100 * paramResults(i,j).FP / paramResults(i,j).TotalNoiseOpportunities;
                         else
                              warning('Visualizer: FP field missing in performance results for params (%d, %d).', i, j);
                         end
                     else
                          % Set to NaN if no noise tests were run for this param set
                          fpr_matrix(i,j) = NaN;
                          % fprintf('Visualizer: No noise opportunities for params (%d, %d), FPR set to NaN.\n', i, j);
                     end
                 end
             end

             figHandle = figure;
             param1Name = obj.MSCAnalyzer.analyzedParam1Name;
             param2Name = obj.MSCAnalyzer.analyzedParam2Name;
             param1Vals = obj.MSCAnalyzer.analyzedParam1Values;
             param2Vals = obj.MSCAnalyzer.analyzedParam2Values;

             if (nParam1 > 1 && nParam2 == 1) || (nParam1 == 1 && nParam2 > 1) || strcmp(p.PlotType, 'line')
                 % --- Line Plot (if one parameter varied or forced) ---
                 if nParam1 > 1 && nParam2 == 1
                     plotData = fpr_matrix(:, 1); % FPR values for varying param1
                     plotParams = param1Vals;
                     xlabelStr = sprintf('%s', param1Name);
                     titleStr = sprintf('FPR vs %s (%s = %.1f)', param1Name, param2Name, param2Vals(1));
                 elseif nParam1 == 1 && nParam2 > 1
                      plotData = fpr_matrix(1, :); % FPR values for varying param2
                      plotParams = param2Vals;
                      xlabelStr = sprintf('%s', param2Name);
                      titleStr = sprintf('FPR vs %s (%s = %.1f)', param2Name, param1Name, param1Vals(1));
                 else % Both > 1, but line plot requested (plot slices)
                     warning('Visualizer: Line plot requested for 2 varying parameters. Plotting FPR vs %s for each %s.', param1Name, param2Name);
                     hold on;
                     colors = lines(nParam2); % Get different colors
                     for j=1:nParam2
                        plot(param1Vals, fpr_matrix(:,j), '-o', 'Color', colors(j,:), 'DisplayName', sprintf('%s = %.2f', param2Name, param2Vals(j)));
                     end
                     hold off;
                      xlabelStr = sprintf('%s', param1Name);
                      titleStr = sprintf('FPR vs %s (Slices for %s)', param1Name, param2Name);
                      legend show;
                      plotData = []; % Skip single plot command below
                 end

                 if ~isempty(plotData)
                     plot(plotParams, plotData, '-o');
                     xlabel(xlabelStr);
                 end
                 ylabel('False Positive Rate (%)');
                 grid on;
                 title(titleStr);

             elseif nParam1 > 1 && nParam2 > 1 && strcmp(p.PlotType, 'image')
                 % --- Image Plot (if both parameters varied) ---
                 imagesc(param2Vals, param1Vals, fpr_matrix); % x=param2, y=param1
                 axis xy; % Set origin to lower-left
                 colorbar;
                 xlabel(sprintf('%s', param2Name));
                 ylabel(sprintf('%s', param1Name));
                 title(sprintf('False Positive Rate (%%) vs Parameters (%s Strategy)', obj.SequentialTester.Strategy.Name));
                 clim([0, max(10, max(fpr_matrix(:)))]); % Set color limits, ensure at least 0-10 range
             else
                  warning('Visualizer: Cannot determine plot type for FPR. nParam1=%d, nParam2=%d, PlotType=%s', nParam1, nParam2, p.PlotType);
                  close(figHandle); figHandle = []; % Close empty figure
             end
        end % plotFalsePositiveRate

        function figHandle = plotParetoFront(obj, p)
            % Plots a Pareto front based on calculated performance metrics.
            % Requires AverageStoppingTime to be calculated and stored.
            arguments
                obj
                p.XAxisMetric char {mustBeMember(p.XAxisMetric, {'FPR', 'FNR'})} = 'FPR'
                p.YAxisMetric char {mustBeMember(p.YAxisMetric, {'AverageStoppingTime'})} = 'AverageStoppingTime'
                p.YAxisUnits {mustBeMember(p.YAxisUnits, {'stages', 'seconds'})} = 'stages' % Units for stopping time
            end

            figHandle = []; % Initialize

            % --- Check if required data exists ---
             if isempty(obj.SequentialTester) || isempty(obj.SequentialTester.PerformanceMetrics) || ...
               ~isfield(obj.SequentialTester.PerformanceMetrics, 'PerParameter')
                 warning('Visualizer: Performance metrics per parameter set not calculated. Cannot plot Pareto front.'); return;
             end
             % Check for Y-axis metric (Stopping Time)
              if strcmp(p.YAxisMetric, 'AverageStoppingTime')
                  if ~isfield(obj.SequentialTester.PerformanceMetrics, 'AverageStoppingTime')
                     warning('Visualizer: Average stopping time not calculated. Run tester.calculateAverageStoppingTime first.'); return;
                  else
                      yDataMatrix = obj.SequentialTester.PerformanceMetrics.AverageStoppingTime;
                       yLabelStr = sprintf('Average Stopping Time (%s)', p.YAxisUnits);
                  end
             else
                 error('Visualizer: Y-axis metric "%s" not currently supported for Pareto plot.', p.YAxisMetric);
             end


             % --- Extract X-axis Metric Data ---
             paramResults = obj.SequentialTester.PerformanceMetrics.PerParameter;
             [nParam1, nParam2] = size(paramResults);
             xDataMatrix = zeros(nParam1, nParam2) * NaN; % Initialize with NaN

             for i = 1:nParam1
                 for j = 1:nParam2
                     res = paramResults(i,j);
                     if strcmp(p.XAxisMetric, 'FPR')
                         if isfield(res, 'TotalNoiseOpportunities') && res.TotalNoiseOpportunities > 0 && isfield(res, 'FP')
                             xDataMatrix(i,j) = 100 * res.FP / res.TotalNoiseOpportunities;
                             xLabelStr = 'False Positive Rate (%)';
                         end
                     elseif strcmp(p.XAxisMetric, 'FNR')
                          if isfield(res, 'TotalSignalOpportunities') && res.TotalSignalOpportunities > 0 && isfield(res, 'FN')
                             xDataMatrix(i,j) = 100 * res.FN / res.TotalSignalOpportunities;
                              xLabelStr = 'False Negative Rate (%)';
                          end
                     end % Add other X metrics if needed
                 end
             end

            % --- Prepare data points for Pareto calculation ---
            % Combine X and Y data into a list of points [x, y]
            % We want to minimize both X and Y (FPR/FNR and Stopping Time)
            validMask = ~isnan(xDataMatrix) & ~isnan(yDataMatrix);
            points = [xDataMatrix(validMask), yDataMatrix(validMask)];
             if isempty(points)
                warning('Visualizer: No valid data points found for Pareto analysis.'); return;
             end

            % Get corresponding parameter indices for labeling
             [rowIdx, colIdx] = find(validMask);
             paramIndices = [rowIdx, colIdx]; % Store indices of valid points

            % --- Find Pareto Front ---
            % Pareto front seeks to minimize objectives. Lower FPR/FNR and lower StopTime are better.
             [paretoPoints, paretoIdx] = obj.findParetoFront(points);
             if isempty(paretoPoints)
                 warning('Visualizer: No Pareto optimal points found.'); return;
             end

             % --- Plot ---
             figHandle = figure; hold on; grid on; box on;
             title(sprintf('Pareto Front: %s vs. %s (%s Strategy)', yLabelStr, xLabelStr, obj.SequentialTester.Strategy.Name));
             xlabel(xLabelStr);
             ylabel(yLabelStr);

             % Plot all valid points
             plot(points(:,1), points(:,2), '.k', 'MarkerSize', 8, 'DisplayName', 'All Parameter Sets');

             % Plot Pareto optimal points
             plot(paretoPoints(:,1), paretoPoints(:,2), 'o-r', 'MarkerSize', 8, 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'Pareto Front');

              % Add labels for Pareto points using parameter values
             param1Name = obj.MSCAnalyzer.analyzedParam1Name;
             param1Vals = obj.MSCAnalyzer.analyzedParam1Values;
             param2Name = obj.MSCAnalyzer.analyzedParam2Name;
             param2Vals = obj.MSCAnalyzer.analyzedParam2Values;

             for k = 1:size(paretoPoints, 1)
                 originalIdx = paretoIdx(k); % Index within the 'points' array
                 pIdx = paramIndices(originalIdx, :); % Get row (p1) and col (p2) index
                 p1Val = param1Vals(pIdx(1));
                 p2Val = param2Vals(pIdx(2));
                 text(paretoPoints(k,1), paretoPoints(k,2) * 0.99, ... % Position label slightly below point
                      sprintf('\\{%s=%.1f, %s=%.1f\\}', param1Name, p1Val, param2Name, p2Val), ...
                      'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8);
             end

             legend show;
             hold off;
        end % plotParetoFront

    end % methods

    methods (Static, Access = private)
        function [ pFront, idxs] = findParetoFront( p )
            % Filters a set of points P according to Pareto dominance for MINIMIZATION.
            % Points that are dominated (non-strictly) by any other point are removed.
            % Inputs:
            % - P    : N-by-D matrix of points (objectives to be minimized).
            % Outputs:
            % - pFront: Pareto-filtered points.
            % - idxs  : Indices of the non-dominated points in the original P.

            [N, D] = size(p);
            isDominated = false(N, 1); % Track dominated points

            for i = 1:N
                if isDominated(i); continue; end % Skip if already marked as dominated

                % Check if point i dominates any other point j
                % Dominates if p(i,d) <= p(j,d) for all d, and p(i,d) < p(j,d) for at least one d
                 dominates_any = false(N, 1);
                 for j = (i+1):N % Only need to check points after i
                     if isDominated(j); continue; end

                     if all(p(i,:) <= p(j,:)) && any(p(i,:) < p(j,:))
                         isDominated(j) = true; % Point j is dominated by i
                     elseif all(p(j,:) <= p(i,:)) && any(p(j,:) < p(i,:))
                          isDominated(i) = true; % Point i is dominated by j
                          break; % No need to check further for point i
                     end
                 end
            end

            pFront = p(~isDominated, :);
            originalIndices = (1:N)';
            idxs = originalIndices(~isDominated);
        end % findParetoFront (Static, Private)
    end % Static Private methods

end % classdef