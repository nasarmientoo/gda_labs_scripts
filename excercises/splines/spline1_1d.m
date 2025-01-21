% Linear Spline Interpolation with Multiple Grid Points (t)

% Step 1: Observations and their Coordinates
points = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6'};
t = [1.8, 2.2, 3.6, 3.8, 6.0, 6.4];
q0 = [9.2, 10.2, 9.3, 8.8, 8.4, 7.9];

disp('Step 1: Observations and Coordinates:');
disp(table(points', t', q0', 'VariableNames', {'Point', 't', 'q0'}));

% Step 2: Reference Grid Points (Multiple Centers)
grid_points = [2,4,6];
disp('Step 2: Grid Points (Spline Centers):');
disp(grid_points');

% Step 3: Compute s(t) for Each Observation and Grid Point
% Compute the spline weights s(t) for each observation relative to each grid point
s1 = @(t) max(0, 1 - abs(t)); % Spline function for linear interpolation
s_t = zeros(length(t), length(grid_points));
for i = 1:length(grid_points)
    s_t(:, i) = arrayfun(@(ti) s1(ti - grid_points(i)), t);
end

disp('Step 3: Spline Weights (s(t)) for Each Grid Point:');
disp(array2table(s_t, 'VariableNames', strcat('Grid_', string(grid_points))));

% Step 4: Design Matrix (A)
% The design matrix for linear splines incorporates all grid points
A = s_t;

disp('Step 4: Design Matrix (A):');
disp(A);

% Step 5: Observations Vector
b = q0'; % Observations vector
disp('Step 5: Observations Vector (b):');
disp(b);

% Step 6: Normal Matrix (N)
N = A' * A;
disp('Step 6: Normal Matrix (N = A^T * A):');
disp(N);

% Step 7: Regularization
% Regularization for multiple grid points
lambda = 1; % Regularization parameter
R = eye(length(grid_points)); % Regularization matrix (identity matrix)
N_reg = N + lambda * R;

disp('Step 7: Regularized Normal Matrix (N_reg):');
disp(N_reg);

% Step 8: Right-Hand Side
rhs = A' * b;
disp('Step 8: Right-Hand Side (A^T * b):');
disp(rhs);

% Step 9: Solve for Coefficients
a = N_reg \ rhs;
disp('Step 9: Spline Coefficients (a):');
disp(table(grid_points', a, 'VariableNames', {'GridPoint', 'Coefficient'}));

% Step 10: Interpolation of New Points
new_points = [0.3, 0.7, 1.0]; % New t-values for interpolation
interpolated_values = zeros(size(new_points));
for i = 1:length(new_points)
    s_new = arrayfun(@(gp) s1(new_points(i) - gp), grid_points);
    interpolated_values(i) = s_new * a;
end

disp('Step 10: Interpolated values for new points:');
disp(table(new_points', interpolated_values', 'VariableNames', {'t', 'InterpolatedValue'}));
