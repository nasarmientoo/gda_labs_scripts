% Spline Interpolation with Least Squares Adjustment

% Step 1: Observations and Time Points
t = [0.9, 1.2, 1.6, 1.7, 2.2, 2.6, 3.4, 3.6]'; % Time points
y0 = [1.26, 1.70, -0.38, -0.47, 1.12, 0.84, 0.99, -0.18]'; % Observations
disp('Step 1: Observations and time points:');
disp('t = '); disp(t);
disp('y0 = '); disp(y0);

% Step 2: Design Matrix
A = [1 0 0; 
     1 0 0; 
     0 1 0; 
     0 1 0; 
     0 1 0; 
     0 0 1; 
     0 0 1;  
     0 0 0];
disp('Step 2: Design matrix (A):');
disp(A);

% Step 3: Normal Matrix
N = A' * A;
disp('Step 3: Normal matrix (N = A^T * A):');
disp(N);

% Step 4: AT * Y0
ATY0 = A' * y0;
disp('Step 4: Product of A^T and Y0 (A^T * Y0):');
disp(ATY0);

% Step 5: Solve for Coefficients
a = inv(N) * ATY0;
disp('Step 5: Coefficients (a = N^-1 * A^T * Y0):');
disp(a);

% Step 6: Residuals
Y_est = A * a;
v = y0 - Y_est;
disp('Step 6: Estimated values (Y_est = A * a):');
disp(Y_est);
disp('Residuals (v = Y0 - A * a):');
disp(v);

% Step 7: Error Variance
n = length(y0); % Number of observations
m = length(a);  % Number of coefficients
sigma0_squared = (v' * v) / (n - m);
sigma0 = sqrt(sigma0_squared);
disp('Step 7: Error variance (sigma0^2 = v^T * v / (n - m)):');
disp(sigma0_squared);
disp('Standard deviation (sigma0):');
disp(sigma0);

% Step 8: Covariance Matrix of Coefficients
C_aa = sigma0_squared * inv(N);
disp('Step 8: Covariance matrix of coefficients (C_aa = sigma0^2 * N^-1):');
disp(C_aa);

% Step 9: Regularization Matrix Construction
disp('Step 9: Regularization matrix (R):');
% The regularization matrix R is derived to enforce smoothness:
% Regularization term: (a_B - a_A)^2 + (a_C - a_B)^2
% Expand and identify contributions:
% R = [  1 -1  0
%        -1  2 -1
%         0 -1  1]
R = [ 1 -1  0;   % Penalty between a_A and a_B
     -1  2 -1;   % Penalty for a_B (connected to both a_A and a_C)
      0 -1  1];  % Penalty between a_B and a_C
disp(R);

% Step 10: Regularized Normal Matrix
N_reg = N + R;
disp('Step 10: Regularized normal matrix (N_reg = N + R):');
disp(N_reg);

% Step 11: Coefficients with Regularization
a_reg = inv(N_reg) * ATY0;
disp('Step 11: Coefficients with regularization (a_reg = N_reg^-1 * A^T * Y0):');
disp(a_reg);

% Step 12: Standard Deviations with Regularization
C_aa_reg = sigma0_squared * inv(N_reg);
sigmaA_reg = sqrt(C_aa_reg(1,1));
sigmaB_reg = sqrt(C_aa_reg(2,2));
sigmaC_reg = sqrt(C_aa_reg(3,3));
disp('Step 12: Standard deviations with regularization:');
disp(['sigmaA_reg = ', num2str(sigmaA_reg)]);
disp(['sigmaB_reg = ', num2str(sigmaB_reg)]);
disp(['sigmaC_reg = ', num2str(sigmaC_reg)]);
