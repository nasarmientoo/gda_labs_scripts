% Practice 1: Least Square interpolation and Gram Schmidt matrix
% orthogonalize 
t = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]';
y0 = [30.4833360027761, 31.9139548258792, 35.0769860573871, 37.3950695229859, ...
      39.9445580803979, 36.9973881232153, 34.6728070449349, 37.93716765258, ...
      47.9497510285442, 47.6384003466922, 49.3976169260511, 44.9932619723902, ...
      52.1392146785302, 50.5095285440181, 50.49573755734, 54.1266939043491, ...
      63.2258909252702, 57.306798540023, 60.0700648208676, 63.7663273700701, ...
      63.6679189134027, 68.5825859319869, 74.1297989258579, 76.3047227897367, ...
      78.6487592765861]';

% Least square solution 1
A = [ones(length(t), 1), t];
N = A'*A
Ninverse = inv(N)
x = Ninverse*A'*y0
y1 = A*x
v = y0 - y1

% Least square solution 2
A2 = [ones(length(t), 1), t, t.^2];
N = A2'*A2
Ninverse = inv(N)
x2 = Ninverse*A2'*y0
y2 = A2*x2
v2 = y0 - y2

% Least square solution 3
A3 = [ones(length(t), 1), t, t.^2, t.^3];
N = A3'*A3
Ninverse = inv(N)
x3 = Ninverse*A3'*y0
y3 = A3*x3
v3 = y0 - y3

%Comparing the error average
mean_v = mean(v);
mean_v2 = mean(v2);
mean_v3 = mean(v3);

% Plot y0 and y1 against t
plot(t, y0, '-o', 'DisplayName', 'y0'); % Plot y0 with markers
hold on; % Keep the plot for adding y1
plot(t, y1, '-s', 'DisplayName', 'y1'); % Plot y1 with different markers
hold on; % Keep the plot for adding y1
plot(t, y2, '-+', 'DisplayName', 'y2'); % Plot y2 with different markers
hold on; % Keep the plot for adding y1
plot(t, y3, '-*', 'DisplayName', 'y3'); % Plot y3 with different markers
hold off;

% Add labels and legend
xlabel('t');
ylabel('y');
title('Plot of y0, y1, y2, and y3 based on t');
legend show; % Display the legend
grid on;     % Optional: Add grid for better readability

% Get the y-axis limits to place the text dynamically
ylimits = ylim;

% Adding mean values dynamically, placing them at 90% of the y-axis range
text(5, ylimits(1) + 0.9*(ylimits(2)-ylimits(1)), ['Mean v: ', num2str(mean_v)], 'FontSize', 10);
text(10, ylimits(1) + 0.9*(ylimits(2)-ylimits(1)), ['Mean v2: ', num2str(mean_v2)], 'FontSize', 10);
text(15, ylimits(1) + 0.9*(ylimits(2)-ylimits(1)), ['Mean v3: ', num2str(mean_v3)], 'FontSize', 10);



