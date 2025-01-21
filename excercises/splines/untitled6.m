% Inputs: Time observations and observed values
t = [1.8, 2.2, 3.6, 3.8, 6.0, 6.4];
y = [9.2, 10.2, 9.3, 8.8, 8.4, 7.9];
spline_centers = [2, 4, 6];

% Step 1: Construct the design matrix for linear splines
A = zeros(length(t), length(spline_centers));
for i = 1:length(spline_centers)
    A(:, i) = max(0, 1 - abs(t - spline_centers(i)) / 2);
end
disp(A)

% Step 2: Solve for spline coefficients
coeffs = (A' * A) \ (A' * y');

% Step 3: Compute residuals and standard deviation
residuals = y - (A * coeffs)';
sigma2_res = sum(residuals.^2) / (length(y) - length(coeffs));
sigma_res = sqrt(sigma2_res);

% Step 4: Interpolate at tP = 3
B = max(0, 1 - abs(3 - spline_centers) / 2);
y_est_tP = B * coeffs;
sigma2_tP = B * (inv(A' * A)) * B' * sigma2_res;