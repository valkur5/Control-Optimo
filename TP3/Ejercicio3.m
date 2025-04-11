%%TP3, fourier de tren de pulsos
clear all, close all, clc;
% Definir parámetros
T = 1;          % Período (segundos)
tau = 0.5;      % Ancho del pulso (segundos)
A = 1;          % Amplitud
fs = 1000;      % Frecuencia de muestreo (Hz)
t_max = 5;      % Tiempo total (segundos)

% Crear vector de tiempo
t = 0:1/fs:t_max;
N = length(t);  % Número de muestras

% Generar el tren de pulsos rectangulares
y = zeros(size(t));  % Inicializar la señal
for i = 1:length(t)
    % Ajustar mod(t, T) para que el pulso esté centrado
    t_mod = mod(t(i) + T/2, T) - T/2;  % Desplazar para centrar el pulso
    if abs(t_mod) <= tau/2
        y(i) = A;
    else
        y(i) = 0;
    end
end

% Calcular la transformada de Fourier
Y = fft(y);  % Transformada rápida de Fourier
f = (0:N-1)*(fs/N);  % Vector de frecuencias
Y_magnitude = abs(Y)/N;  % Magnitud normalizada
Y_phase = angle(Y);      % Fase

% Ajustar el espectro para frecuencias positivas
f = f(1:floor(N/2));
Y_magnitude = Y_magnitude(1:floor(N/2));
Y_phase = Y_phase(1:floor(N/2));
Y_magnitude(2:end) = 2 * Y_magnitude(2:end);  % Duplicar magnitudes (excepto DC) por simetría

% Frecuencia fundamental
f0 = 1/T;  % Frecuencia fundamental (Hz)

% Reconstruir la señal con las componentes cosenoidales
num_harmonics = 10;  % Número de armónicos a considerar (ajustable)
y_reconstructed = zeros(size(t));  % Señal reconstruida

% Componente DC (frecuencia 0)
y_reconstructed = y_reconstructed + Y_magnitude(1) * cos(2 * pi * 0 * t + Y_phase(1));

% Sumar los armónicos
for k = 1:num_harmonics
    idx = find(abs(f - k*f0) == min(abs(f - k*f0)), 1);  % Encontrar índice de la frecuencia k*f
    if ~isempty(idx)
        amplitude = Y_magnitude(idx);
        phase = Y_phase(idx);
        y_reconstructed = y_reconstructed + amplitude * cos(2 * pi * (k*f0) * t + phase);
    end
end

% Graficar la señal original y la reconstruida superpuestas
figure 1;
plot(t, y, 'LineWidth', 1, 'DisplayName', 'Señal Cuadrada Original');
hold on; grid on;
plot(t, y_reconstructed, 'LineWidth', 1, 'DisplayName', sprintf('Reconstrucción (%d armónicos)', num_harmonics));
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Cuadrada y Reconstrucción con Componentes Cosenoidales');
legend('show');

