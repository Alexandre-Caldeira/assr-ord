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
            % Check if data exists at the specified indices
            if isempty(mscData) && (idx1 > size(obj.MSCAnalyzer.mscValues,1) || idx2 > size(obj.MSCAnalyzer.mscValues,2) || p1Idx > size(obj.MSCAnalyzer.mscValues,3) || p2Idx > size(obj.MSCAnalyzer.mscValues,4))
                 warning('Visualizer: Indices (Idx1:%d, Idx2:%d, P1:%d, P2:%d) exceed bounds of calculated MSC data. Cannot plot.', idx1, idx2, p1Idx, p2Idx); return;
            elseif isempty(mscData)
                 warning('Visualizer: MSC data is empty for specified indices (Idx1:%d, Idx2:%d, P1:%d, P2:%d). Cannot plot.', idx1, idx2, p1Idx, p2Idx); return;
            end

            decisions = obj.SequentialTester.TestDecisions{idx1, idx2, p1Idx, p2Idx};
            stopping = obj.SequentialTester.StoppingStages{idx1, idx2, p1Idx, p2Idx};
            strategyParams = obj.SequentialTester.StrategyParamsArray{idx1, idx2, p1Idx, p2Idx};
            epochs = obj.MSCAnalyzer.getEpochs(idx1, idx2, p1Idx, p2Idx);
            M_vals = obj.MSCAnalyzer.getWindowsPerStage(idx1, idx2, p1Idx, p2Idx);

            if isempty(decisions) || isempty(stopping) || isempty(epochs) || isempty(M_vals) || isempty(strategyParams)
                 warning('Visualizer: Decisions, stopping stages, epochs, M values, or strategyParams missing for the specified indices (Idx1:%d, Idx2:%d, P1:%d, P2:%d). Cannot plot.', idx1, idx2, p1Idx, p2Idx); return;
            end

            [nBins, K, nChan] = size(mscData);
            if isempty(nChan) || nChan == 0 % Handle case where mscData might be 2D
                if ndims(mscData) == 2 && channelIdx == 1
                   nChan = 1; % Assume single channel if data is 2D and channel 1 requested
                else
                     warning('Visualizer: MSC data has unexpected dimensions or 0 channels. Cannot plot.'); return;
                end
            end
            if channelIdx < 1 || channelIdx > nChan
                error('Visualizer: Selected channel index %d is out of bounds [1, %d].', channelIdx, nChan);
            end
            if K == 0; warning('Visualizer: MSC data has 0 stages. Cannot plot.'); return; end

            % --- Recalculate Thresholds ---
            analysisInfo.K = K;
            analysisInfo.M = M_vals;
            if numel(epochs)>1; analysisInfo.MM = epochs(2:end) - 1; else; analysisInfo.MM=[]; end % Required for ETS
            try
                thresholds = obj.SequentialTester.Strategy.calculateThresholds(strategyParams, analysisInfo);
            catch ME
                 warning('Visualizer: Could not recalculate thresholds for plotting. Error: %s', ME.message); return;
            end

            % --- Prepare Data for Plotting ---
            mscSequenceChannel = squeeze(mscData(:, :, channelIdx));
            if K==1; mscSequenceChannel = mscSequenceChannel(:); end % Ensure column vector if K=1
            cumulativeMSC = cumsum(mscSequenceChannel, 2);

            % Identify signal/noise bins
            sigFreqBins = []; noiseFreqBins = [];
            if strcmp(obj.DataManager.Mode, 'sim')
                 fs = obj.DataManager.Fs; nfft = obj.DataManager.getNFFT(); maxBin = obj.DataManager.getNumBins();
                 if nfft > 0
                     freqResolution = fs / nfft;
                     sigFreqs_Hz = obj.DataManager.SimSignalFrequencies; noiseFreqs_Hz = obj.DataManager.SimNoiseFrequencies;
                     sigFreqBins = round(sigFreqs_Hz / freqResolution) + 1; noiseFreqBins = round(noiseFreqs_Hz / freqResolution) + 1;
                     sigFreqBins = unique(sigFreqBins(sigFreqBins >= 1 & sigFreqBins <= maxBin));
                     noiseFreqBins = unique(noiseFreqBins(noiseFreqBins >= 1 & noiseFreqBins <= maxBin));
                     noiseFreqBins = setdiff(noiseFreqBins, sigFreqBins);
                 else
                      warning('Visualizer: NFFT is 0, cannot map frequencies to bins.');
                 end
            else; sigFreqBins = 1:nBins; end % Default for exp mode

            % --- Create Plot ---
            figHandle = figure; hold on; grid on; box on;
            param1Name = obj.MSCAnalyzer.analyzedParam1Name; param1ValStr = '?';
            if ~isempty(obj.MSCAnalyzer.analyzedParam1Values) && p1Idx <= numel(obj.MSCAnalyzer.analyzedParam1Values); param1ValStr = mat2str(obj.MSCAnalyzer.analyzedParam1Values(p1Idx)); end
            param2Name = obj.MSCAnalyzer.analyzedParam2Name; param2ValStr = '?';
             if ~isempty(obj.MSCAnalyzer.analyzedParam2Values) && p2Idx <= numel(obj.MSCAnalyzer.analyzedParam2Values); param2ValStr = mat2str(obj.MSCAnalyzer.analyzedParam2Values(p2Idx)); end

            plotTitle = sprintf('MSC vs. Thresholds (%s=%s, %s=%s)\nIdx1=%d, Idx2=%d, Chan=%d, Strategy=%s', ...
                               param1Name, param1ValStr, param2Name, param2ValStr, ...
                               idx1, idx2, channelIdx, obj.SequentialTester.Strategy.Name);
            title(plotTitle); xlabel('Stage (k)');

             if strcmp(obj.SequentialTester.Strategy.Name,'CGST_Beta'); ylabel('Cumulative MSC'); plotData = cumulativeMSC;
             else; ylabel('Per-Stage MSC'); plotData = mscSequenceChannel; end % ETS uses per-stage

            hPlots = []; legendEntries = {};

            if ~isempty(sigFreqBins) % Plot Signal Bins
                 validSigBins = sigFreqBins(sigFreqBins <= nBins);
                 if ~isempty(validSigBins)
                    hPlots(end+1) = plot(1:K, plotData(validSigBins, :)', '-b');
                    legendEntries{end+1} = 'Signal Freq Bins';
                 end
            end
            if ~isempty(noiseFreqBins) % Plot Noise Bins
                 validNoiseBins = noiseFreqBins(noiseFreqBins <= nBins);
                  if ~isempty(validNoiseBins)
                    hPlots(end+1) = plot(1:K, plotData(validNoiseBins, :)', ':r');
                    legendEntries{end+1} = 'Noise Freq Bins';
                  end
             end

            % Plot Thresholds
            if isfield(thresholds, 'efficacy') && ~isempty(thresholds.efficacy)
                hPlots(end+1) = plot(1:K, thresholds.efficacy, '--g', 'LineWidth', 2);
                legendEntries{end+1} = 'Efficacy Threshold';
            end
             if isfield(thresholds, 'detection') && ~isempty(thresholds.detection)
                hPlots(end+1) = plot(1:K, thresholds.detection, '--g', 'LineWidth', 2);
                legendEntries{end+1} = 'Detection Threshold';
                if isfield(thresholds,'NDC') && thresholds.NDC > 1; legendEntries{end} = sprintf('Detection Threshold (NDC=%d)', thresholds.NDC); end
            end
            if obj.SequentialTester.Strategy.RequiresFutilityThreshold && isfield(thresholds, 'futility') && ~isempty(thresholds.futility)
                hPlots(end+1) = plot(1:K, thresholds.futility, '--m', 'LineWidth', 2);
                legendEntries{end+1} = 'Futility Threshold';
            end

            if ~isempty(hPlots) % Add legend
                [~, uniqueIdx] = unique(hPlots); legend(hPlots(uniqueIdx), legendEntries(uniqueIdx));
            end

            % Adjust Y-axis limits
            allYData = plotData(:); % Use the data actually plotted (cumulative or per-stage)
            if isfield(thresholds, 'efficacy') && ~isempty(thresholds.efficacy); allYData = [allYData; thresholds.efficacy(:)]; end
            if isfield(thresholds, 'detection') && ~isempty(thresholds.detection); allYData = [allYData; thresholds.detection(:)]; end
            if isfield(thresholds, 'futility') && ~isempty(thresholds.futility); allYData = [allYData; thresholds.futility(:)]; end
            validYData = allYData(~isnan(allYData) & ~isinf(allYData));
            if ~isempty(validYData); minY = min(validYData); maxY = max(validYData); ylim([min(0, minY - 0.1*abs(minY)), max(maxY + 0.1*abs(maxY), 0.1)]);
            else; ylim([0, 1]); end % Default if no valid data
            xlim([0.5, K+0.5]);

            hold off;
        end % plotMSCandThresholds

        function figHandle = plotFalsePositiveRate(obj, p)
            % Plots False Positive Rate (FPR) across the parameter sets tested.
            arguments
                obj
                p.PlotType {mustBeMember(p.PlotType, {'line', 'image'})} = 'image'
            end
            figHandle = [];
             if isempty(obj.SequentialTester) || isempty(obj.SequentialTester.PerformanceMetrics) || ~isfield(obj.SequentialTester.PerformanceMetrics, 'PerParameter')
                 warning('Visualizer: Performance metrics not calculated. Cannot plot FPR. Run calculatePerformance first.'); return;
             end
             paramResults = obj.SequentialTester.PerformanceMetrics.PerParameter;
             if isempty(paramResults); warning('Visualizer: No PerParameter performance results found.'); return; end
             [nParam1, nParam2] = size(paramResults);

             fpr_matrix = nan(nParam1, nParam2); % Initialize with NaN
             for i = 1:nParam1; for j = 1:nParam2
                 if isfield(paramResults(i,j), 'TotalNoiseOpportunities') && paramResults(i,j).TotalNoiseOpportunities > 0 && isfield(paramResults(i,j), 'FP')
                     fpr_matrix(i,j) = 100 * paramResults(i,j).FP / paramResults(i,j).TotalNoiseOpportunities;
                 end; end; end

             figHandle = figure;
             param1Name = obj.MSCAnalyzer.analyzedParam1Name; param1Vals = obj.MSCAnalyzer.analyzedParam1Values;
             param2Name = obj.MSCAnalyzer.analyzedParam2Name; param2Vals = obj.MSCAnalyzer.analyzedParam2Values;

             if (nParam1 > 1 && nParam2 == 1) || (nParam1 == 1 && nParam2 > 1) || strcmp(p.PlotType, 'line')
                 if nParam1 > 1 && nParam2 == 1; plotData = fpr_matrix(:, 1); plotParams = param1Vals; xLabelStr = param1Name; titleStr = sprintf('FPR vs %s (%s=%s)', param1Name, param2Name, mat2str(param2Vals(1)));
                 elseif nParam1 == 1 && nParam2 > 1; plotData = fpr_matrix(1, :); plotParams = param2Vals; xLabelStr = param2Name; titleStr = sprintf('FPR vs %s (%s=%s)', param2Name, param1Name, mat2str(param1Vals(1)));
                 else % Both > 1, line plot slices
                     warning('Visualizer: Plotting FPR slices vs %s for each %s.', param1Name, param2Name);
                     hold on; colors = lines(nParam2);
                     for j=1:nParam2; plot(param1Vals, fpr_matrix(:,j), '-o', 'Color', colors(j,:), 'DisplayName', sprintf('%s = %s', param2Name, mat2str(param2Vals(j)))); end
                     hold off; xLabelStr = param1Name; titleStr = sprintf('FPR vs %s (Slices for %s)', param1Name, param2Name); legend show; plotData = [];
                 end
                 if ~isempty(plotData); plot(plotParams, plotData, '-o'); xlabel(xLabelStr); end
                 ylabel('False Positive Rate (%)'); grid on; title(titleStr);
             elseif nParam1 > 1 && nParam2 > 1 && strcmp(p.PlotType, 'image')
                 imagesc(param2Vals, param1Vals, fpr_matrix); axis xy; colorbar;
                 xlabel(param2Name); ylabel(param1Name);
                 title(sprintf('False Positive Rate (%%) vs Parameters (%s Strategy)', obj.SequentialTester.Strategy.Name));
                 clim([0, max([10, max(fpr_matrix(:),[],'omitnan')])]); % Ensure reasonable color limits
             else; warning('Visualizer: Cannot determine plot type for FPR.'); close(figHandle); figHandle = []; end
        end % plotFalsePositiveRate

        function figHandle = plotParetoFront(obj, p)
            % Plots a Pareto front based on calculated performance metrics.
            arguments
                obj
                p.XAxisMetric char {mustBeMember(p.XAxisMetric, {'FPR', 'FNR'})} = 'FPR'
                p.YAxisMetric char {mustBeMember(p.YAxisMetric, {'AverageStoppingTime'})} = 'AverageStoppingTime'
                p.YAxisUnits {mustBeMember(p.YAxisUnits, {'stages', 'seconds'})} = 'stages'
            end
            figHandle = [];
             if isempty(obj.SequentialTester) || isempty(obj.SequentialTester.PerformanceMetrics) || ~isfield(obj.SequentialTester.PerformanceMetrics, 'PerParameter')
                 warning('Visualizer: PerParameter performance metrics not calculated.'); return; end
             if strcmp(p.YAxisMetric, 'AverageStoppingTime') && ~isfield(obj.SequentialTester.PerformanceMetrics, 'AverageStoppingTime')
                 warning('Visualizer: Average stopping time not calculated.'); return;
             elseif strcmp(p.YAxisMetric, 'AverageStoppingTime'); yDataMatrix = obj.SequentialTester.PerformanceMetrics.AverageStoppingTime; yLabelStr = sprintf('Average Stopping Time (%s)', p.YAxisUnits);
             else; error('Visualizer: Y-axis metric "%s" not supported.', p.YAxisMetric); end

             paramResults = obj.SequentialTester.PerformanceMetrics.PerParameter;
             [nParam1, nParam2] = size(paramResults);
             xDataMatrix = nan(nParam1, nParam2);

             for i = 1:nParam1; for j = 1:nParam2; res = paramResults(i,j);
                 if strcmp(p.XAxisMetric, 'FPR') && isfield(res,'TotalNoiseOpportunities') && res.TotalNoiseOpportunities > 0 && isfield(res,'FP')
                     xDataMatrix(i,j) = 100 * res.FP / res.TotalNoiseOpportunities; xLabelStr = 'False Positive Rate (%)';
                 elseif strcmp(p.XAxisMetric, 'FNR') && isfield(res,'TotalSignalOpportunities') && res.TotalSignalOpportunities > 0 && isfield(res,'FN')
                     xDataMatrix(i,j) = 100 * res.FN / res.TotalSignalOpportunities; xLabelStr = 'False Negative Rate (%)';
                 end; end; end

            validMask = ~isnan(xDataMatrix) & ~isnan(yDataMatrix);
            points = [xDataMatrix(validMask), yDataMatrix(validMask)];
            if isempty(points); warning('Visualizer: No valid points for Pareto analysis.'); return; end
            [rowIdx, colIdx] = find(validMask); paramIndices = [rowIdx, colIdx];

             [paretoPoints, paretoIdx] = obj.findParetoFront(points); % Assumes minimization
             if isempty(paretoPoints); warning('Visualizer: No Pareto optimal points found.'); return; end

             figHandle = figure; hold on; grid on; box on;
             title(sprintf('Pareto Front: %s vs. %s (%s Strategy)', yLabelStr, xLabelStr, obj.SequentialTester.Strategy.Name));
             xlabel(xLabelStr); ylabel(yLabelStr);
             plot(points(:,1), points(:,2), '.k', 'MarkerSize', 8, 'DisplayName', 'All Parameter Sets');
             plot(paretoPoints(:,1), paretoPoints(:,2), 'o-r', 'MarkerSize', 7, 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', 'Pareto Front');

             param1Name = obj.MSCAnalyzer.analyzedParam1Name; param1Vals = obj.MSCAnalyzer.analyzedParam1Values;
             param2Name = obj.MSCAnalyzer.analyzedParam2Name; param2Vals = obj.MSCAnalyzer.analyzedParam2Values;
             for k = 1:size(paretoPoints, 1)
                 originalIdx = paretoIdx(k); pIdx = paramIndices(originalIdx, :);
                 p1ValStr = mat2str(param1Vals(pIdx(1)),3); % Use mat2str for flexible formatting
                 p2ValStr = mat2str(param2Vals(pIdx(2)),3);
                 text(paretoPoints(k,1), paretoPoints(k,2), ...
                      sprintf(' \\{%s=%s, %s=%s\\}', param1Name(1), p1ValStr, param2Name(1), p2ValStr), ... % Abbreviate param names
                      'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'FontSize', 7);
             end
             legend show Location best; hold off;
        end % plotParetoFront

        function figHandle = plotConfusionMatrix(obj, p1Idx, p2Idx)
            % Plots a confusion matrix for a specific parameter set index.
            % Requires simulation mode and calculated performance metrics.
            % Requires Deep Learning Toolbox for confusionchart.
             figHandle = [];
             if isempty(obj.SequentialTester) || isempty(obj.SequentialTester.PerformanceMetrics) || ~isfield(obj.SequentialTester.PerformanceMetrics, 'PerParameter')
                 warning('Visualizer: PerParameter performance metrics not calculated. Cannot plot confusion matrix.'); return;
             end
              if ~strcmp(obj.DataManager.Mode, 'sim')
                  warning('Visualizer: Confusion matrix requires simulation mode ground truth.'); return;
              end
               if ~license('test', 'Deep_Learning_Toolbox')
                  warning('Visualizer: Deep Learning Toolbox not found. Cannot plot confusion matrix using confusionchart.'); return;
               end

             paramResults = obj.SequentialTester.PerformanceMetrics.PerParameter;
             [nParam1, nParam2] = size(paramResults);
             if p1Idx < 1 || p1Idx > nParam1 || p2Idx < 1 || p2Idx > nParam2
                 warning('Visualizer: Parameter indices (%d, %d) out of bounds [%d x %d].', p1Idx, p2Idx, nParam1, nParam2); return;
             end

             res = paramResults(p1Idx, p2Idx);
             TP = res.TP; FP = res.FP; TN = res.TN; FN = res.FN;
             TotalP = res.TotalSignalOpportunities; % Ground Truth Positive (Signal)
             TotalN = res.TotalNoiseOpportunities;  % Ground Truth Negative (Noise)

             % Check if counts are valid
             if isempty(TP) || isempty(FP) || isempty(TN) || isempty(FN) || TotalP == 0 || TotalN == 0
                 warning('Visualizer: Performance counts invalid or missing for parameter set (%d, %d). Cannot plot confusion matrix.', p1Idx, p2Idx); return;
             end

              % Create the 2x2 confusion matrix
             % Rows: True Class (Signal, Noise)
             % Columns: Predicted Class (Signal(Detect=1), Noise(Reject=0 or -1))
             % CM = [TP FN; FP TN]
             cmData = [TP, FN; FP, TN];

             figHandle = figure;
             cm = confusionchart(cmData, {'Signal', 'Noise'}, {'Predicted Signal', 'Predicted Noise'});

             param1Name = obj.MSCAnalyzer.analyzedParam1Name; param1ValStr = '?';
             if ~isempty(obj.MSCAnalyzer.analyzedParam1Values) && p1Idx <= numel(obj.MSCAnalyzer.analyzedParam1Values); param1ValStr = mat2str(obj.MSCAnalyzer.analyzedParam1Values(p1Idx),3); end
             param2Name = obj.MSCAnalyzer.analyzedParam2Name; param2ValStr = '?';
             if ~isempty(obj.MSCAnalyzer.analyzedParam2Values) && p2Idx <= numel(obj.MSCAnalyzer.analyzedParam2Values); param2ValStr = mat2str(obj.MSCAnalyzer.analyzedParam2Values(p2Idx),3); end

             cm.Title = sprintf('Confusion Matrix (%s=%s, %s=%s) - %s', ...
                                param1Name, param1ValStr, param2Name, param2ValStr, obj.SequentialTester.Strategy.Name);
             % Optional: Customize appearance (e.g., normalization)
             % cm.Normalization = 'row-normalized'; % Show percentages (TPR, FNR; FPR, TNR)
             % cm.Normalization = 'column-normalized'; % Show precision, false discovery rate, etc.
             sortClasses(cm, {'Signal', 'Noise'}); % Ensure consistent order
        end % plotConfusionMatrix


        function figHandle = plotOverallPerformanceBars(obj)
            % Plots a bar chart of the overall TPR, FPR, TNR, FNR.
            % Requires simulation mode and calculated performance metrics.
             figHandle = [];
             if isempty(obj.SequentialTester) || isempty(obj.SequentialTester.PerformanceMetrics) || ~isfield(obj.SequentialTester.PerformanceMetrics, 'Overall')
                 warning('Visualizer: Overall performance metrics not calculated. Cannot plot bars.'); return;
             end
              if ~strcmp(obj.DataManager.Mode, 'sim')
                  warning('Visualizer: Overall performance bars require simulation mode ground truth.'); return;
              end

              overallPerf = obj.SequentialTester.PerformanceMetrics.Overall;
              metrics = {'TPR', 'FNR', 'TNR', 'FPR'}; % Order for plotting
              values = nan(1, 4);
              if isfield(overallPerf, 'TPR'); values(1) = overallPerf.TPR; end
              if isfield(overallPerf, 'FNR'); values(2) = overallPerf.FNR; end
              if isfield(overallPerf, 'TNR'); values(3) = overallPerf.TNR; end
              if isfield(overallPerf, 'FPR'); values(4) = overallPerf.FPR; end

              if all(isnan(values))
                   warning('Visualizer: All overall performance metric values are NaN. Cannot plot bars.'); return;
              end

              figHandle = figure;
              bar(values);
              set(gca, 'xticklabel', metrics);
              ylabel('Rate (%)');
              ylim([0 100]);
              grid on;
              title(sprintf('Overall Performance Metrics (%s Strategy)', obj.SequentialTester.Strategy.Name));

              % Add text labels to bars
               textHandles = text(1:numel(values), values, sprintfc('%.1f%%', values), ...
                                  'HorizontalAlignment','center', 'VerticalAlignment','bottom');
               % Adjust text position slightly if needed, e.g., for short bars
                for i = 1:numel(textHandles)
                    if values(i) < 5 % Adjust label position for very small bars
                         textHandles(i).VerticalAlignment = 'top';
                         textHandles(i).Position(2) = textHandles(i).Position(2) - 1; % Shift down slightly
                    end
                end

        end % plotOverallPerformanceBars


    end % methods

    methods (Static, Access = private)
        function [ pFront, idxs] = findParetoFront( p )
            % Filters a set of points P according to Pareto dominance for MINIMIZATION.
            [N, ~] = size(p);
            isDominated = false(N, 1);
            for i = 1:N
                if isDominated(i); continue; end
                 for j = (i+1):N
                     if isDominated(j); continue; end
                     if all(p(i,:) <= p(j,:)) && any(p(i,:) < p(j,:)); isDominated(j) = true; % i dominates j
                     elseif all(p(j,:) <= p(i,:)) && any(p(j,:) < p(i,:)); isDominated(i) = true; break; end % j dominates i
                 end
            end
            pFront = p(~isDominated, :);
            originalIndices = (1:N)';
            idxs = originalIndices(~isDominated);
        end % findParetoFront (Static, Private)
    end % Static Private methods

end % classdef