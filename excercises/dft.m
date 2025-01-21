% Script unificado para Análisis de Fourier y Predicciones
% Define los datos directamente en el script.

% Paso 1: Observaciones y Puntos de Tiempo
disp('Paso 1: Observaciones y Puntos de Tiempo');

% Ejemplo de datos: Modifica estas líneas según el ejercicio

y0 = [13.3, 17.4, 10.2, 2.4, -3.7, -0.2, 7.5, -1.5, -12.1, -2.6 ]'; 
n = length(y0); % Número de observaciones
t = (0:n-1)'; % 
if length(t) ~= length(y0)
    error('Los puntos de tiempo (t) y las observaciones (y0) deben tener la misma longitud.');
end


disp('Datos de entrada:');
disp(table(t, y0, 'VariableNames', {'Tiempo', 'Observaciones'}));

% Verificar regularidad de los puntos de tiempo
if any(diff(t) ~= mean(diff(t)))
    disp('Los puntos de tiempo son irregularmente espaciados.');
else
    disp('Los puntos de tiempo son regularmente espaciados.');
end

% Paso 2: Calcular Coeficientes Armónicos
disp('Paso 2: Cálculo de Coeficientes Armónicos usando DFT.');
a0 = (1/n) * sum(y0); % Coeficiente de orden 0
disp(['Fórmula: a0 = (1/n) * sum(y0), Resultado: a0 = ', num2str(a0)]);

% Determinar índice máximo armónico basado en n (par o impar)
% if mod(n, 2) == 0
%     k_max = n / 2; % Para n par
% else
%     k_max = (n - 1)/ 2; % Para n impar
% end
% disp('Frecuencia máxima')
% disp(k_max / n)
k_max=2
% Inicializar matrices para coeficientes
ak = zeros(k_max, 1);
bk = zeros(k_max, 1);

% Calcular ak y bk para todas las frecuencias
for k = 1:k_max
    ak(k) = (1/n) * sum(y0 .* cos(2 * pi * t * k / n));
    bk(k) = (1/n) * sum(y0 .* sin(2 * pi * t * k / n));
end


disp('Coeficientes armónicos calculados:');
disp(table((1:k_max)', ak, bk, 'VariableNames', {'Frecuencia', 'ak', 'bk'}));

% Paso 3: Desviación Estándar del Error
sigma0_squared = (sum(y0.^2) - n * a0^2 - 2 * n * sum(ak.^2 + bk.^2)) / (n - 2 * k_max - 1);
sigma0 = sqrt(sigma0_squared);
disp(['Fórmula: sigma0^2 = (sum(y0^2) - n*a0^2 - 2*n*sum(ak^2 + bk^2)) / (n - 2*k_max - 1)']);
disp(['Resultado: sigma0 = ', num2str(sigma0)]);

% Paso 4: Frecuencia Dominante
amplitudes = sqrt(ak.^2 + bk.^2);
[~, dominant_idx] = max(amplitudes);
period = n / dominant_idx; % Periodo del armónico dominante
disp(['La frecuencia dominante tiene un periodo de: ', num2str(period), ' unidades de tiempo.']);

% Paso 5: Reconstruir Señal Usando el Modelo Armónico
disp('Paso 5: Reconstrucción de la señal usando el modelo armónico.');
t_interp = t; % Puntos de tiempo para la interpolación
y_hat = a0 + zeros(size(t_interp));

for k = 1:k_max
    y_hat = y_hat + 2 * (ak(k) * cos(2 * pi * t_interp * k / n) + bk(k) * sin(2 * pi * t_interp * k / n));
end

disp('Señal reconstruida:');
disp(table(t_interp, y_hat, 'VariableNames', {'Tiempo', 'Valor Reconstruido'}));

% Paso 6: Predicción para un Tiempo Futuro
future_t = 10; % Modifica este valor para predecir en otro tiempo
y_future = a0;

for k = 1:k_max
    y_future = y_future + 2 * (ak(k) * cos(2 * pi * future_t * k / n) + bk(k) * sin(2 * pi * future_t * k / n));
end

disp(['Valor predicho en t = ', num2str(future_t), ': ', num2str(y_future)]);

% Error de Predicción
prediction_variance = sigma0_squared * (1 + sum(1 ./ amplitudes.^2));
prediction_error = sqrt(prediction_variance);
disp(['Error de predicción en t = ', num2str(future_t), ': ', num2str(prediction_error)]);
disp(prediction_variance)
