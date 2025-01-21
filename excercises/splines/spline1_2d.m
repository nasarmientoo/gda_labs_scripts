% Bilinear Spline Interpolation

% Step 1: Observations and their Coordinates
points = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'};
x = [0.2, 0.8, 0.4, 1.2, 1.4, 0.1, 0.6, 0.4]; % x-coordinates
y = [0.0, 1.2, 1.0, 0.4, 0.9, 0.4, 0.8, 1.4]; % y-coordinates
q0 = [8.7, 9.8, 7.8, 8.3, 9.7, 8.9, 9.4,7.9]; % Observations
disp('Step 1: Observations and Coordinates:');
table(points', x', y', q0', 'VariableNames', {'Point', 'x', 'y', 'q0'})


% Step 2: Grid Points (Spline Centers)
grid_points = {'PA', 'PB', 'PC', 'PD'};
x_grid = [0, 0, 1, 1]; % x-coordinates of grid points
y_grid = [0, 1, 0, 1]; % y-coordinates of grid points
disp('Step 2: Grid Points (Spline Centers):');
table(grid_points', x_grid', y_grid', 'VariableNames', {'GridPoint', 'x', 'y'})

% Step 3: Matrix of s(x) and s(y) Values for Each Point
% Compute the individual s(x) and s(y) values for each observation point and grid point
s1 = @(t) max(0, 1 - abs(t)); % Spline function
s_x = zeros(length(q0), length(grid_points));
s_y = zeros(length(q0), length(grid_points));

for i = 1:length(q0)
    for j = 1:length(grid_points)
        s_x(i, j) = s1((x(i) - x_grid(j)) / 2);
        s_y(i, j) = s1((y(i) - y_grid(j)) / 2);
    end
end

disp('Step 3: Matrix of s(x) values for each observation and grid point:');
disp(s_x);
disp('Step 3: Matrix of s(y) values for each observation and grid point:');
disp(s_y);

% Step 4: Design Matrix (A)
% Multiply s(x) and s(y) to compute the design matrix A
A = s_x .* s_y;

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

%%
% Step 7: Regularization Matrix
% Adding smoothness constraints in the y-direction
D = [0, 1, -1, 0;
    1, 0, 0, -1;];
R = D'*D; % Regularization matrix
lambda = 1; % Regularization parameter
N_reg = N + lambda * R;

disp('Step 7: Regularized Normal Matrix (N_reg):');
disp(R)
disp(N_reg);

% Step 8: Right-Hand Side
rhs = A' * b;
disp('Step 8: Right-Hand Side (A^T * b):');
disp(rhs);

% Step 9: Solve for Coefficients
a = N_reg \ rhs;
disp('Step 9: Spline Coefficients (a):');
disp(a);


