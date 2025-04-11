#TP3, graficas de las autocorrelaciones}

%%Señal triangular (Autocorrelación de onda cuadrada)
T = 4; A = 1; fs = 1000; t_max = 10;
t = 0:1/fs:t_max;
y_triang = A * (1 - 2 * abs(mod(t, T) - T/2) / T);

figure 1;
plot(t, y_triang); xlabel('\tau'); ylabel('Amplitud'); title('\Phi_y_y'); grid on;

%%Señal sinc (Autocorrelacion de ruido blanco filtrado)
y_sync=sinc(t);
figure 2;
plot(t,y_sync);xlabel('\tau'); ylabel('Amplitud'); title('\Phi_x_x'); grid on;

%%Señal constante (Autocorrelación de una señal constante)
y_const=ones(size(t));
figure 3;
plot(t,y_const);xlabel('\tau'); ylabel('Amplitud'); title('\Phi_z_z'); grid on;

%%Salida suma (Autocorrelación w de la suma de todas las anteriores)
y_sal=y_triang.+y_sync.+y_const;
figure 4;
plot(t,y_sal);xlabel('\tau'); ylabel('Amplitud'); title('\Phi_w_w'); grid on;

%%Señal senoidal
figure 5;
plot(t,A*sin(t-pi/4)); grid on; xlabel("tiempo [s]");ylabel("Amplitud"); title("5*sin(t-\pi/4)")

%%Autocorrelacion de la señal senoidal
plot(t,A*[-sin(t-pi/4)]+A); xlabel('\tau'); ylabel('Amplitud'); title('\Phi_x_x'); grid on;
