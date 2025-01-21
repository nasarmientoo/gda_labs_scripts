function weightedAverageCalculations()
    % Given data points (x, y, z, q)
    points = [
        1, 0, 0, 6.2;  % P1
        0, 1, 1, 7.4;  % P2
        0, 1, 0, 6.8;  % P3
        1, 1, 1, 7.6   % P4
    ];
    
    % Estimate value at (0, 0, 0) first
    estimateAtOrigin(points);

    % Perform Leave-One-Out Validation
    [rmsError, estimatedValues] = weightedAverageLeaveOneOut(points);

    % Display the RMS error
    disp(['RMS Error (LOO): ' num2str(rmsError)]);

end

function [rmsError, estimatedValues] = weightedAverageLeaveOneOut(points)

    n = size(points, 1); % Number of known points
    estimatedValues = zeros(n, 1); % Pre-allocate to hold all estimated values


    % Leave-One-Out Loop
    for i = 1:n
        % Extract estimation point coordinates
        xp = points(i, 1);
        yp = points(i, 2);
        zp = points(i, 3);

        % Create the data points for the current LOO step, with the i-th point
        % excluded.
        looPoints = points;
        looPoints(i,:) = [];

        n_loo = size(looPoints, 1); % Number of known points (excluding 1)
        weights = zeros(n_loo, 1); % Pre-allocate weights
        distances = zeros(n_loo, 1); % Pre-allocate distances

       % Calculate the Euclidean distances and weights for each point, using
       % the current estimation coordinates and the other data points in looPoints
        for j = 1:n_loo
            x = looPoints(j, 1);
            y = looPoints(j, 2);
            z = looPoints(j, 3);
            distance = sqrt((xp - x)^2 + (yp - y)^2 + (zp - z)^2);
            distances(j) = distance;  % Store the distance
            weights(j) = 1 / (2 * distance^2 + 1);
        end

        % Calculate the weighted average with data excluding the i-th point
        weightedSum = sum(weights .* looPoints(:, 4));
        totalWeight = sum(weights);
        estimatedValues(i) = weightedSum / totalWeight;

        % Display intermediate values for each LOO step
        fprintf('LOO step %d:\n', i);
        disp('Distances to other points:');
        for k = 1:n_loo
           fprintf('Point %d distance: %.4f\n',k, distances(k));
        end
        disp('Weights:');
        for k= 1:n_loo
            fprintf('Point %d weight: %.4f\n', k, weights(k));
        end

         fprintf('Estimated value (excluding point %d): %.4f\n\n', i, estimatedValues(i));
    end

    % Calculate the RMS error
    squaredErrors = (estimatedValues - points(:, 4)).^2;
    meanSquaredError = mean(squaredErrors);
    rmsError = sqrt(meanSquaredError);

end


function estimateAtOrigin(points)
    % Coordinates of the point to estimate
    xp = 0;
    yp = 0;
    zp = 0;

    n = size(points, 1);  % Number of known points
    weights = zeros(n, 1);  % Pre-allocate weights
    distances = zeros(n, 1); % Pre-allocate distances


    % Calculate the Euclidean distances and weights for each point
    for j = 1:n
        x = points(j, 1);
        y = points(j, 2);
        z = points(j, 3);
        distance = sqrt((xp - x)^2 + (yp - y)^2 + (zp - z)^2);
        distances(j) = distance;  % Store the distance
         weights(j) = 1 / (2 * distance^2 + 1);
    end

    % Calculate the weighted average
    weightedSum = sum(weights .* points(:, 4));
    totalWeight = sum(weights);
    estimatedValue = weightedSum / totalWeight;

    % Display results
    fprintf('Estimation at P(0, 0, 0):\n');
    disp('Distances to other points:');
     for k = 1:n
        fprintf('Point %d distance: %.4f\n',k, distances(k));
     end
     disp('Weights:');
     for k = 1:n
          fprintf('Point %d weight: %.4f\n', k, weights(k));
     end
    fprintf('Estimated value at (0, 0, 0): %.4f\n\n', estimatedValue);

end