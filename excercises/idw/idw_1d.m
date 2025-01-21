% Define the observation times and values
t_obs = [0.9, 1.2, 1.6, 1.7, 2.2, 2.6, 3.4, 3.6];
y_obs = [1.26, 1.70, -0.38, -0.47, 1.12, 0.84, 0.99, -0.18];

% Define the grid nodes
t_grid = [1, 2, 3];

% Initialize arrays to store distances, weights, and estimates
tau = zeros(length(t_grid), length(t_obs));
weights = zeros(length(t_grid), length(t_obs));
y_est_grid = zeros(1, length(t_grid)); % Estimates for grid nodes
y_est_loo = zeros(1, length(t_obs));  % Estimates using leave-one-out method

% Calculate weights and estimates for grid nodes
for i = 1:length(t_grid)
    numerator = 0;
    denominator = 0;
    for j = 1:length(t_obs)
        % Calculate the time distance
        tau(i, j) = abs(t_grid(i) - t_obs(j));
        
        % Calculate the weight based on the given function
        if tau(i, j) <= 0.5
            weights(i, j) = 1 - 2 * tau(i, j);
        else
            weights(i, j) = 0;
        end
        
        % Accumulate numerator and denominator for grid estimates
        numerator = numerator + weights(i, j) * y_obs(j);
        denominator = denominator + weights(i, j);
    end
    % Calculate the estimate y_est for the current grid node
    if denominator ~= 0
        y_est_grid(i) = numerator / denominator;
    else
        y_est_grid(i) = NaN;
    end
end

% Perform leave-one-out method for all 8 points
y_est_loo = zeros(1, length(t_obs));
weights_loo = zeros(length(t_obs), length(t_obs)); % Store LOO weights

for i = 1:length(t_obs)
    numerator_loo = 0;
    denominator_loo = 0;
    for j = 1:length(t_obs)
        if j ~= i % Exclude the current observation
            tau_loo = abs(t_obs(i) - t_obs(j));
            if tau_loo <= 0.5
                weights_loo(i, j) = 1 - 2 * tau_loo;
            else
                weights_loo(i, j) = 0;
            end
            numerator_loo = numerator_loo + weights_loo(i, j) * y_obs(j);
            denominator_loo = denominator_loo + weights_loo(i, j);
        end
    end
    % Calculate the leave-one-out estimate
    if denominator_loo ~= 0
        y_est_loo(i) = numerator_loo / denominator_loo;
    else
        y_est_loo(i) = NaN;
    end
end

% Calculate RMS error for grid nodes
rms_grid = sqrt(mean((y_est_grid - y_obs(1:length(t_grid))).^2, 'omitnan'));

% Calculate RMS error for leave-one-out estimates
rms_loo = sqrt(mean((y_est_loo - y_obs).^2, 'omitnan'));

% Display results with formulas
disp('-------------------------- Formulas Used --------------------------');
disp('1. Weighting Function:');
disp('   w(τ) = 1 - 2τ, for 0 ≤ τ ≤ 0.5');
disp('   w(τ) = 0, for τ > 0.5');
disp(' ');
disp('2. Estimate Formula (for each node):');
disp('   y_est = Σ(w(i) * y(i)) / Σ(w(i))');
disp(' ');
disp('3. RMS Error Formula:');
disp('   RMS = sqrt(mean((y_est - y_obs)^2))');
disp('-------------------------------------------------------------------');

disp('Distances (tau) for all grid points:');
disp(array2table(tau, 'VariableNames', strcat("t_", string(1:length(t_obs))), ...
                      'RowNames', strcat("Grid_", string(1:length(t_grid)))));

disp('Weights for all grid points:');
disp(array2table(weights, 'VariableNames', strcat("t_", string(1:length(t_obs))), ...
                        'RowNames', strcat("Grid_", string(1:length(t_grid)))));

disp('Estimated values at grid nodes (y_est_grid):');
disp(y_est_grid);

disp('Leave-one-out weights for all points:');
disp(array2table(weights_loo, 'VariableNames', strcat("t_", string(1:length(t_obs))), ...
                             'RowNames', strcat("LOO_Point_", string(1:length(t_obs)))));

disp('Estimated values with leave-one-out method (y_est_loo):');
disp(y_est_loo);

disp('RMS error for grid nodes:');
disp(['RMS_grid = ', num2str(rms_grid)]);

disp('RMS error for leave-one-out method:');
disp(['RMS_LOO = ', num2str(rms_loo)]);
