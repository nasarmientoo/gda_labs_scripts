function weightedAverageCalculations2D()
    % Given data points (x, y, q)
    points = [
        1, 0, 10.2;  
        2, 2, 9.2;   
        0, 1, 9.8;  
        0, 2, 7.8    
    ];

    % Estimate value at (0, 0) and perform Leave-One-Out Validation
    estimateValueAndLOO(points, [0, 0]);
end

function estimateValueAndLOO(points, targetPoint)
    % Estimate value at the target point
    estimatedValue = calculateWeightedAverage(points, targetPoint);
    fprintf('Estimated value at (%0.1f, %0.1f): %.4f\n\n', targetPoint(1), targetPoint(2), estimatedValue);

    % Perform Leave-One-Out Validation
    n = size(points, 1);
    estimatedValues = zeros(n, 1);

    for i = 1:n
        looPoints = points;
        looPoints(i, :) = [];  % Exclude the i-th point
        estimatedValues(i) = calculateWeightedAverage(looPoints, points(i, 1:2));
    end

    % Calculate and display RMS error
    rmsError = sqrt(mean((estimatedValues - points(:, 3)).^2));
    disp('Weighted Average: sqrt(sum(estimatedValues * values))^2 / n');
    disp(['RMS Error (LOO): ' num2str(rmsError)]);
end

function estimatedValue = calculateWeightedAverage(points, targetPoint)
    % Extract target coordinates
    xp = targetPoint(1);
    yp = targetPoint(2);

    % Calculate distances and weights
    distances = sqrt((points(:, 1) - xp).^2 + (points(:, 2) - yp).^2);
    weights = 1 ./ (distances.^2 + 1);

    % Compute weighted average
    weightedSum = sum(weights .* points(:, 3));
    totalWeight = sum(weights);
    estimatedValue = weightedSum / totalWeight;

    % Display intermediate results
    disp('Distances and Weights:');
    for i = 1:length(distances)
        fprintf('Point %d: Distance = %.4f, Weight = %.4f\n', i, distances(i), weights(i));
    end
end
