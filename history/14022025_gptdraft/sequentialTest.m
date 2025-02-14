function [decisions, stageMetrics] = sequentialTest(MSCvalues, aThresholds, gThresholds, K) 
% sequentialTest applies the decision process using the accumulated MSC values.
    % decisions = zeros(1, K); 
    stageMetrics = cumsum(MSCvalues,2); 
    decisions = zeros(size(stageMetrics));

    for k = 1:K 
        decisions(stageMetrics(:,k,:) >= aThresholds(k)) = 1;
        decisions(stageMetrics(:,k,:) < gThresholds(k)) = -1;

        % if stageMetrics(:,k,:) >= aThresholds(k) 
        % 
        %     decisions(k) = 1; 
        % 
        %     % Detection: reject H0 
        %     decisions(k+1:end) = NaN; 
        % 
        %     % Early stopping 
        %     break; 
        % 
        % elseif stageMetrics(k) <= gThresholds(k) 
        % 
        %     decisions(k) = -1; 
        % 
        %     % Futility: accept H0 (stop test) 
        %     decisions(k+1:end) = NaN; 
        % 
        %     break; 
        % 
        % else 
        % 
        %     decisions(k) = 0; 
        %     % Continue to next stage 
        % 
        % end 
    end 
end

%

% onde stageMetrics(:,k,:) >= aThresholds(k), stageMetrics(:,k,:) = 1