% M - Observations
% N - Parameters
% System with N equations with M unknowns
% y - Noiseless Data KNOWN
% t - Puntos conocidos en el eje t
% q - Función base

function a = exactInterpolationSystem(y, q, t)

n = length(t);

% Construir la matriz Q con las funciones base
Q = zeros(n, n);
for i = 1:n
    for j = 1:n
        Q(i, j) = q{j}(t(i));
    end
end

% Resolver el sistema de ecuaciones: Q * a = y
a = Q \ y;
end

% Ejemplo de uso:
% Nodos conocidos y valores:
t = [1, 2, 3];
y = [4, 5, 6];

% Definir funciones base q1, q2, q3 como ejemplo:
q = {@(t) 1, @(t) t, @(t) t.^2};

% Calcular los coeficientes de interpolación:
a = exactInterpolationSystem(y, q, t);

% Mostrar los coeficientes:
fprintf('Los coeficientes de interpolación son:\n');
for i = 1:length(a)
    fprintf('a%d = %.2f\n', i, a(i));
end

% Reconstrucción de y utilizando los coeficientes calculados:
Q = zeros(length(t), length(t));
for i = 1:length(t)
    for j = 1:length(t)
        Q(i, j) = q{j}(t(i));
    end
end

reconstructed_y = Q * a;

% Mostrar la reconstrucción:
fprintf('\nReconstrucción de y:\n');
for i = 1:length(reconstructed_y)
    fprintf('y_reconstructed(%d) = %.2f\n', i, reconstructed_y(i));
end