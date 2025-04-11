%%TP3, Ejercicio 4
clear all, close all, clc;

K=1;T=1;
Fw=tf([K],[T 1])
figure 1; bode(Fw);

t=0:0.01:10;

figure 2; hold on; grid on;
for i=1:5
  plot(t,K/i*exp(-t/i),"LineWidth",2);
end
xlabel("tiempo [s]"); ylabel("amplitud");
title("Funcion de correlacion, K=1")
legend("T=1","T=2","T=3","T=4","T=5")
