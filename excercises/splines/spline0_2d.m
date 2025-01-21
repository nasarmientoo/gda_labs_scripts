% Zero-Order Spline Interpolation in 2D

% Step 1: Observations and their Coordinates
points = {'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8'};
x = [0.2, 0.8, 0.4, 1.2, 1.4, 0.1, 0.6, 0.4]; % x-coordinates
y = [0.0, 1.2, 1.0, 0.4, 0.9, 0.4, 0.8, 1.4]; % y-coordinates
q0 = [8.7, 9.8, 7.8, 8.3, 9.7, 8.9, 9.4, 7.9]; % Observations
disp('Step 1: Observations and Coordinates:');
table(points', x', y', q0', 'VariableNames', {'Point', 'x', 'y', 'q0'})

% Step 2: Grid Points (Spline Centers)
grid_points = {'PA', 'PB', 'PC', 'PD'};
x_grid = [0, 0, 1, 1]; % x-coordinates of grid points
y_grid = [0, 1, 0, 1]; % y-coordinates of grid points
disp('Step 2: Grid Points (Spline Centers):');
table(grid_points', x_grid', y_grid', 'VariableNames', {'GridPoint', 'x', 'y'})

% Step 3: Assign Each Observation to the Closest Grid Point
num_obs = length(q0);
num_grid = length(grid_points);
A = zeros(num_obs, num_grid);

for i = 1:num_obs
    % Compute distances to each grid point
    distances = sqrt((x(i) - x_grid).^2 + (y(i) - y_grid).^2);
    % Find the closest grid point
    [~, closest] = min(distances);
    % Assign 1 to the closest grid point column
    A(i, closest) = 1;
end

disp('Step 3: Design Matrix (A):');
disp(A);

% Step 4: Observations Vector
b = q0'; % Observations vector
disp('Step 4: Observations Vector (b):');
disp(b);

% Step 5: Normal Matrix (N)
% If the noise standard deviation is correlated then   N = inv(A'* inv(Q) *
% A) and (a = inv(N) * A^T inv(Q)* Y0)  if not N = A'A and (a = N^-1 * A^T * Y0) 
% Q is the covariance matrix, with variances along the diagonal.
Q = [1, 0, 0, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0, 0;
     0, 0, 0, 0, 4, 0, 0, 0;
     0, 0, 0, 0, 0, 4, 0, 0;
     0, 0, 0, 0, 0, 0, 4, 0;
     0, 0, 0, 0, 0, 0, 0, 4];
N = (A' / Q) * A;
disp('Step 5: Normal Matrix (N = A^T * A):');
disp(N);

% Step 6: Solve for Coefficients
%rhs = A' * b;
%a = N \ rhs;  Solve the normal equations
a = N \ (A' / Q) * b;
disp('Step 6: Spline Coefficients (a):');
disp(a);

% Step 7: Display the Interpolation Function
disp('Interpolation Function:');
for i = 1:num_grid
    fprintf('  f(x, y) = %.3f, for points closest to Grid Point (%d, %d)\n', ...
        a(i), x_grid(i), y_grid(i));
end

R = [1, -1, 0, 0;
     0, 1, 0, -1;
     0, 0, -1, 1;
    -1, 0, 1, 0;];
D = N + 0.1 * (R' * R)

disp(D)

% Step 8: Visualization
figure;
hold on;
grid on;
title('Zero-Order Spline Interpolation in 2D');
xlabel('x');
ylabel('y');

% Plot the grid points
scatter(x_grid, y_grid, 100, 'k', 'filled', 'DisplayName', 'Grid Points');

% Plot the observations
scatter(x, y, 50, 'r', 'filled', 'DisplayName', 'Observations');

% Annotate the regions (constant values)
for i = 1:num_grid
    text(x_grid(i), y_grid(i), sprintf('%.3f', a(i)), 'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'right', 'Color', 'blue', 'FontSize', 12);
end

legend('Location', 'best');
hold off;
