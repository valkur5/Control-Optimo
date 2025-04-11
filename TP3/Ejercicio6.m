%%TP3, Ejercicio 6
N=10000;
rand('state',0)
x=0;xv=zeros(N,1);
for hh=1:N
    x=((pi+x)^5)-fix((pi+x)^5);
    xv(hh)=x;
end
% hist(xv)
%xv=rand(N,1);
fixx = xcorr(xv,xv);

figure 1;
plot(fixx(N:2*N-1));title('Autocorrelación de x');
xlabel('Tiempo tao')

SxM=fft(fixx(N:2*N-1));

figure 2;
semilogx(20*log10(abs(SxM(1:N/2))),'.k'); grid on;
xlabel('Pulsación en rad por seg')
ylabel("S_{xx}(w)")
title('Módulo de F(w) en dB')
