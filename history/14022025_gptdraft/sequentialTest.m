function [decisions, stageMetrics] = sequentialTest(MSCvalues, params) 
% sequentialTest applies the decision process using the accumulated MSC values.
    % decisions = zeros(1, K); 
    stageMetrics = cumsum(MSCvalues,2); 
    decisions = zeros(size(stageMetrics));

    % % Compute stopping decisions
    % for k = 1:params.K
    %     decisions(stageMetrics(:,k,:) >= params.aThresholds(k)) = 1;
    %     decisions(stageMetrics(:,k,:) < params.gThresholds(k)) = -1;
    % end 

    for channel=1:params.nChannels
        for k = 1:params.K
            for freq = params.testFrequencies

                if ~(k==1) && ~isnan(decisions(channel, k, freq))%(decisions(channel, k, freq))~=0
                    % Not first test and already decided, keep decision
                    decisions(channel, k, freq) = decisions(channel, k-1, freq);

                else
                    if stageMetrics(channel, k, freq) >= params.aThresholds(k)
                        % Enough evidence to detect signal
                        decisions(channel, k, freq) = 1;

                    elseif stageMetrics(channel, k, freq) < params.gThresholds(k)
                        % Enough evidence to stop trial (futile)
                        decisions(channel, k, freq) = -1;

                    elseif k==params.K 
                        % Last trial, and not enough evidence found for
                        % either detection criteria (ERROR!)
                        % decisions(channel, :, freq) = NaN;
                        error('Oops')
                        
                    end

                end

            end
        end
    end


end

%
% onde stageMetrics(:,k,:) >= aThresholds(k), stageMetrics(:,k,:) = 1
%
% onde decisions(:,k+1,:) = decisions(:,1:k,:) if not zero