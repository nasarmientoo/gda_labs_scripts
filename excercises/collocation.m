% Define the observation times and values
t_obs = [1:10];
y_obs = [-3.2, -2.3, -1.3, -2.2, 1.4, 2.4, 3.2, 7.5, 1.0, 2.4];
N = length(t_obs);

% Define the single grid node
t_grid = 2;

% Standard deviation of the noise
sigma_noise = 0.5;

% Covariance function
covariance = @(tau) (tau <= 0.4) .* (1 - abs(tau)/4);


% Initialize Cxy and Cxx
Cxy = zeros(1, length(t_obs)); % Single row for one grid node
Cxx = zeros(length(t_obs), length(t_obs));

% Compute Cxy (covariance between the grid node and observations)
for j = 1:length(t_obs)
    tau = abs(t_grid - t_obs(j)); % Time distance
    Cxy(1, j) = covariance(tau);
end

% Compute Cxx (covariance between observations)
for i = 1:length(t_obs)
    for j = 1:length(t_obs)
        tau = abs(t_obs(i) - t_obs(j)); % Time distance
        Cxx(i, j) = covariance(tau);
    end
end

% Add noise variance to the diagonal of Cxx
Cxx = Cxx + sigma_noise^2 * eye(length(t_obs));

% Compute Yest using the formula Yest = Cxy * inv(Cxx) * Y0
Cxx_inv = inv(Cxx);
y_est = Cxy * Cxx_inv * y_obs';

% Compute the variance of the prediction
Cy0 = covariance(0); % Covariance at zero lag
sigma_est_squared = Cy0 - Cxy * Cxx_inv * Cxy';
sigma_est = sqrt(sigma_est_squared); % Standard deviation


%%
% Predefined tau values for empirical covariance 
tau_values = [0, 1, 2, 3];

% Initialize matrix to store empirical covariances
C_emp_values = zeros(length(tau_values), 1);

% Loop through tau values and compute empirical covariance
for i = 1:length(tau_values)
    tau = tau_values(i);
    % Compute empirical covariance
    if tau == 0
        % Special case for tau = 0
        C_emp = sum(y_obs.^2) / N;
    else
        % General case for tau > 0
        C_emp = sum(y_obs(1:N-tau) .* y_obs(1+tau:N)) / (N - tau);
    end
    
    % Store the result
    C_emp_values(i) = C_emp;
    % Display the result
    disp(['Empirical covariance for tau = ', num2str(tau), ' is: ', num2str(C_emp)]);
end

% Display the design matrix for tau = [1, 2, 3]
A = [1 - 1/4; 1 - 2/4; 1 - 3/4];
disp('Design matrix (A):');
disp(A);

% Compute signal variance (\sigma_y^2) using least squares
C_emp_subset = C_emp_values(2:end); % Exclude tau = 0
AT_A = A' * A;
AT_C = A' * C_emp_subset;
sigma_y2 = AT_C / AT_A;

disp(['Signal variance (\sigma_y^2): ', num2str(sigma_y2)]);

% Compute noise variance (\sigma_\nu^2)
sigma_nu2 = C_emp_values(1) - sigma_y2;
disp(['Noise variance (\sigma_\nu^2): ', num2str(sigma_nu2)]);

% Display the formulas used
disp('Formula for signal variance: \sigma_y^2 = (A^T * A)^{-1} * A^T * C_emp');
disp('Formula for noise variance: \sigma_\nu^2 = C_emp(0) - \sigma_y^2');

%%

% Display formulas
disp('-------------------------- Formulas Used --------------------------');
disp('1. Covariance Function:');
disp('   C_y(τ) = 1 - 2τ, for 0 ≤ τ ≤ 0.5; otherwise, C_y(τ) = 0');
disp(' ');
disp('2. Updated Observation Covariance Matrix:');
disp('   Cxx = Cxx + σ^2 * I');
disp(' ');
disp('3. Predicted Value for t = 2:');
disp('   Yest = Cxy * inv(Cxx) * Y0');
disp(' ');
disp('4. Variance of Prediction for t = 2:');
disp('   σ^2_est = Cy(0) - Cxy * inv(Cxx) * Cxy''');
disp('-------------------------------------------------------------------');

% Display matrices and results
disp('Covariance matrix Cxy (between t = 2 and observations):');
disp(array2table(Cxy, 'VariableNames', strcat("Obs_", string(1:length(t_obs)))));

disp('Covariance matrix Cxx (between observations, with noise variance):');
disp(array2table(Cxx, 'VariableNames', strcat("Obs_", string(1:length(t_obs))), ...
                      'RowNames', strcat("Obs_", string(1:length(t_obs)))));

disp('Inverse of covariance matrix Cxx (Cxx^-1):');
disp(array2table(Cxx_inv, 'VariableNames', strcat("Obs_", string(1:length(t_obs))), ...
                          'RowNames', strcat("Obs_", string(1:length(t_obs)))));

disp('Predicted value (y_est) for t = 2:');
disp(y_est);

disp('Variance of prediction (σ^2_est) for t = 2:');
disp(sigma_est_squared);

disp('Standard deviation (σ_est) for t = 2:');
disp(sigma_est);
