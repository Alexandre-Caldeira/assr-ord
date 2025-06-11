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

        for j = 1:N % Compare with ALL others (can optimize slightly)
           if i == j || isDominated(j)
               continue;
           end

            % Check if j dominates i (assumes minimization)
            if all(data(j, :) <= data(i, :)) && any(data(j, :) < data(i, :))
                isDominated(i) = true;
                break; % Stop checking others if i is dominated
            end
        end
    end

    % Indices of non-dominated points
    frontIndices = find(~isDominated);
    if nargout > 1 % Only calculate points if requested
        frontPoints = data(frontIndices, :);
    end

end