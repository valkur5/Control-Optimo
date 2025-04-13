%%TP3, Ejercicio 7 PRBS9
%pkg load signal;
clear all; close all;

N=500;m=11;
y=zeros(N,1);
x=zeros(1,m+1);x(m1)=1;
for k=1:N
    res=xor(x(3),x(1)); %PRSB11 x^11+x^9+1
    y(k)=2*res-1;
    x_d=circshift(x,[1,1]);
    x_d(1)=res;
    x=x_d;
end
fixx = xcorr(y,y);
SxM=fft(fixx(N:2*N-1));
x

figure 1;
subplot(3,1,1);plot(y,"lineWidth",2); title("PRBS11");ylim([-1.2 1.2]),xlim([0 16]); grid on;
subplot(3,1,2);plot(fixx(N:2*N-1));title('Autocorrelación de x'); xlabel('Tiempo tao'); xlim([0 16]);grid on;
subplot(3,1,3);semilogx(20*log10(abs(SxM(1:N/2))),'.k');xlabel('Pulsación en rad por seg');title('Módulo de F(w) en dB'); grid on;
