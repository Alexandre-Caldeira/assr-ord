% --- FILE START: TestStrategy.m ---
classdef (Abstract) TestStrategy < handle
    % Abstract interface defining the contract for sequential test strategies.
    % Classes implementing this interface can calculate thresholds and
    % determine test outcomes based on a sequence of test statistics (e.g., MSC).

    properties (Abstract)
        % Boolean indicating if this strategy inherently calculates/requires
        % a futility threshold in addition to an efficacy/detection threshold.
        RequiresFutilityThreshold logical

        % Strategy name (e.g., 'CGST_Beta', 'ETS_NDC') used for identification.
        Name char
    end

    methods (Abstract)
        % Calculate or retrieve necessary thresholds for the test.
        % Args:
        %   obj: The strategy object itself.
        %   strategyParams: Struct containing parameters specific to this strategy
        %                   needed for threshold calculation (e.g., .alpha, .NDC, .VC_not).
        %   analysisInfo: Struct containing information about the input data stages,
        %                 needed for context (e.g., .K = number of stages,
        %                 .M = windows per stage (scalar/vector), .MM = window counts per test point).
        % Returns:
        %   thresholds: Struct containing the calculated thresholds (e.g.,
        %               .efficacy, .futility, .detection). Content depends on strategy.
        thresholds = calculateThresholds(obj, strategyParams, analysisInfo);

        % Run the sequential test on a sequence of test statistics (e.g., MSC).
        % Args:
        %   obj: The strategy object itself.
        %   statSequence: The sequence of test statistics (e.g., MSC values) for each stage.
        %                 Expected dimensions: (bins x stages x channels) or similar.
        %   thresholds: The threshold structure returned by calculateThresholds.
        % Returns:
        %   decisionSequence: Matrix indicating the decision at each stage for each bin/channel.
        %                     Values typically: 1 (Detect/Efficacy), -1 (Stop/Futility), 0 (Continue).
        %                     Dimensions usually match statSequence.
        %   stoppingStage: Matrix indicating the stage number at which a decision (1 or -1)
        %                  was first reached for each bin/channel. Value is K if no early stop.
        %                  Dimensions usually: (bins x channels).
        [decisionSequence, stoppingStage] = runTest(obj, statSequence, thresholds);
    end
end