function [frontIndices, frontPoints] = paretoFront(data)
%PARETOFRONT calculates the Pareto front indices and points for a dataset.
%   [frontIndices, frontPoints] = PARETOFRONT(data)
%
%   Input:
%       data: An N-by-M matrix where N is the number of data points and M
%             is the number of objectives. Assumes minimization for all
%             objectives. If maximization is desired for an objective,
%             negate that column before passing it to this function.
%
%   Output:
%       frontIndices: A column vector containing the row indices of the
%                     points belonging to the Pareto front in the original
%                     'data' matrix.
%       frontPoints:  A K-by-M matrix containing the data points that form
%                     the Pareto front (where K is the number of points on
%                     the front).
%
%   Algorithm:
%       A point X dominates point Y if X is strictly better than Y in at
%       least one objective and no worse than Y in all other objectives.
%       The Pareto front consists of all non-dominated points.
%       This implementation uses pairwise comparisons.
%
%   Example:
%       data = [1 5; 2 4; 3 3; 4 2; 5 1; 3 4; 4 3]; % Minimize both objectives
%       [idx, pts] = paretoFront(data);
%       disp('Pareto Front Indices:'); disp(idx);
%       disp('Pareto Front Points:'); disp(pts);
%       figure; plot(data(:,1), data(:,2), 'bo'); hold on;
%       plot(pts(:,1), pts(:,2), 'r*-', 'LineWidth', 1.5);
%       xlabel('Objective 1'); ylabel('Objective 2'); title('Pareto Front Example');
%       legend('All Points', 'Pareto Front'); grid on;

    [N, M] = size(data);
    if N <= 1
        frontIndices = (1:N)'; % If 0 or 1 point, it's on the front
        frontPoints = data;
        return;
    end

    isDominated = false(N, 1); % Initialize: assume no point is dominated

    % Compare each point i with every other point j
    for i = 1:N
        if isDominated(i)
            continue; % Skip if already known to be dominated
        end

        for j = (i+1):N % Compare with subsequent points
            if isDominated(j)
                continue; % Skip if j is already known to be dominated
            end

            % Check for dominance relationship between i and j
            % Assumes minimization!
            i_dominates_j = false;
            j_dominates_i = false;

            % Check if i dominates j
            if all(data(i, :) <= data(j, :)) && any(data(i, :) < data(j, :))
                i_dominates_j = true;
            end

            % Check if j dominates i
            if all(data(j, :) <= data(i, :)) && any(data(j, :) < data(i, :))
                j_dominates_i = true;
            end

            % Update dominated status
            if i_dominates_j
                isDominated(j) = true;
            elseif j_dominates_i
                isDominated(i) = true;
                break; % Stop checking against others if i is dominated
            end
        end
    end

    % Indices of non-dominated points
    frontIndices = find(~isDominated);
    frontPoints = data(frontIndices, :);

end