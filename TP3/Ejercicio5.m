%%TP3, ejercicio 5
pkg load signal;
clear all; close all; clc

fs = 1000; t_max = 10; sigma = 1;
t = 0:1/fs:t_max; N = length(t);

%%Generación de ruido blanco
x = sigma * randn(1, N);

phi_xx=0;tau=0;

[phi_xx, lags] = xcorr(x, 'biased'); tau = lags / fs;
S_xx = fft(phi_xx); f = (-N/2:N/2-1)*(fs/N);
figure;
subplot(2, 1, 1); plot(tau, phi_xx, 'b-', 'LineWidth', 1); xlabel('Desplazamiento \tau (s)'); ylabel('\phi_{xx}(\tau)'); title('Función de Autocorrelación del Ruido Blanco'); grid on; xlim([0 0.5]);
subplot(2, 1, 2);
semilogx(20*log10(abs(S_xx(1:N/2))),'.k');
xlabel('Frecuencia (Hz)'); ylabel('S_{xx}(f)'); title('Autoespectro de Potencia del Ruido Blanco'); grid on; xlim([-fs/2 fs/2]);
