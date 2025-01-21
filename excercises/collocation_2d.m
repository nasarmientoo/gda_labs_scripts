function exercise4()
    % --- Exercise 4: Collocation Prediction and Error Calculation ---

    % --- 1. Define Known Points and Observations ---
    % Coordinates of known points
    P1 = [0, 0];   % Point 1: (x, y)
    P2 = [0, 0.8];   % Point 2: (x, y)
    P3 = [1, 0];   % Point 3: (x, y)
    P4 = [0.8, 1.8];   % Point 3: (x, y)
    points = [P1; P2; P3; P4]; % Matrix with all point locations

    q_observed = [2.8; 3.2; -3.6;-2.4 ]; % Observed values at P1, P2, P3

    % Coordinates of the prediction point
    PA = [0, 0.5];  

    % Noise variance (given)
    sigma_v_squared = 2;

    % --- 2. Define Semi-Variogram and Covariance Functions ---
    % Semi-variogram function
    semiVariogram = @(rho) 1 - exp(-rho.^2);

     % Covariance function
    covariance = @(rho) exp(-rho.^2);

    % --- 3. Function to calculate distance ---
    calculateDistance = @(p1, p2) sqrt(sum((p1 - p2).^2));

    % --- 4. Calculate All Distances ---
    % Distances between known points
     fprintf('\n--- 4. Calculate All Distances --- \n');
     fprintf('   Formula: ρ_ij = √((x_i - x_j)² + (y_i - y_j)²) \n');
    rho12 = calculateDistance(P1, P2);
     fprintf('     ρ12 = distance(P1, P2): %.4f\n', rho12);
    rho13 = calculateDistance(P1, P3);
     fprintf('     ρ13 = distance(P1, P3): %.4f\n', rho13);
    rho23 = calculateDistance(P2, P3);
     fprintf('     ρ23 = distance(P2, P3): %.4f\n', rho23);


    % Distances between known points and prediction point PA
    rho1A = calculateDistance(P1, PA);
     fprintf('     ρ1A = distance(P1, PA): %.4f\n', rho1A);
    rho2A = calculateDistance(P2, PA);
      fprintf('     ρ2A = distance(P2, PA): %.4f\n', rho2A);
    rho3A = calculateDistance(P3, PA);
      fprintf('     ρ3A = distance(P3, PA): %.4f\n', rho3A);

    % --- 5. Construct Covariance Matrix C_XX ---
    % Covariance matrix of known points
    fprintf('\n--- 5. Construct Covariance Matrix C_XX --- \n');
    fprintf('  Formula:\n');
     fprintf('   C_XX =  [C(0) + σv²     C(ρ12)        C(ρ13)   ]\n');
     fprintf('           [ C(ρ21)    C(0) + σv²       C(ρ23)   ]\n');
     fprintf('           [ C(ρ31)     C(ρ32)     C(0) + σv²] \n');

    C_XX = [covariance(0) + sigma_v_squared, covariance(rho12), covariance(rho13);
            covariance(rho12), covariance(0) + sigma_v_squared, covariance(rho23);
            covariance(rho13), covariance(rho23), covariance(0) + sigma_v_squared];
     fprintf('    Where C(ρ) = e^(-ρ²) , C(0) = 1, σv² = %.2f\n', sigma_v_squared);
    fprintf('  C_XX = \n');
    disp(C_XX);

    % --- 6. Construct Covariance Vector C_XA ---
    % Covariance vector between known points and prediction point PA
    fprintf('\n--- 6. Construct Covariance Vector C_XA --- \n');
    fprintf('  Formula:\n');
    fprintf('    C_XA =  [C(ρ1A)]\n');
     fprintf('            [C(ρ2A)]\n');
    fprintf('            [C(ρ3A)]\n');
    fprintf('    Where C(ρ) = e^(-ρ²) \n');
    C_XA = [covariance(rho1A);
            covariance(rho2A);
            covariance(rho3A)];
    fprintf(' C_XA = \n');
    disp(C_XA);

    % --- 7. Calculate Weights (w) ---
    % Calculate the weights using matrix inversion
    fprintf('\n--- 7. Calculate Weights (w) ---\n');
    fprintf('   Formula: w = C_XA'' C_XX ^-1\n');
    w = C_XA' / C_XX; % Equivalent to C_XA' * inv(C_XX)
    fprintf(' w = \n');
    disp(w);


    % --- 8. Predict q at Point PA ---
    % Weighted average of known 'q' values using the weights
    fprintf('\n--- 8. Predict q at Point PA ---\n');
    fprintf('   Formula: q_hat_PA = w * q_observed\n');
    q_hat_PA = w * q_observed;
    fprintf(' Predicted q value (q_hat_PA): %.4f\n', q_hat_PA);


     % --- 9. Calculate Collocation Error Variance ---
    % Error variance of the prediction
      fprintf('\n--- 9. Calculate Collocation Error Variance ---\n');
     fprintf(' Formula: σ²[q_hat(PA)] = C(0) - w * C_XA\n');
     sigma_q_hat_PA_squared = covariance(0) - w * C_XA;
     fprintf('Error Variance (sigma_q_hat_PA_squared): %.4f\n', sigma_q_hat_PA_squared);


    % --- 10. Calculate Error Standard Deviation ---
    % Error standard deviation is the square root of the error variance
     fprintf('\n--- 10. Calculate Error Standard Deviation ---\n');
      fprintf('  Formula: σ[q_hat(PA)] = sqrt(σ²[q_hat(PA)])\n');
    sigma_q_hat_PA = sqrt(sigma_q_hat_PA_squared);
    fprintf(' Error Standard Deviation (sigma_q_hat_PA): %.4f\n', sigma_q_hat_PA);

    % --- 11. Plotting Semi-Variograms ---
   % Create distances for plotting
    rho_values = linspace(0, 3, 100);

    % Calculate semi-variogram values
    semi_variogram_signal_values = semiVariogram(rho_values);
     semi_variogram_obs_values = semiVariogram(rho_values) + sigma_v_squared;

    % Plot the semi-variograms
    figure;
    plot(rho_values, semi_variogram_signal_values, 'b', 'DisplayName', 'Semi-variogram of Signal');
    hold on;
    plot(rho_values, semi_variogram_obs_values, 'r', 'DisplayName', 'Semi-variogram of Observations');
    hold off;
    xlabel('Distance (ρ)');
    ylabel('Semi-variogram (γ)');
    title('Semi-variograms of Signal and Observations');
    legend('show');
    grid on;

end