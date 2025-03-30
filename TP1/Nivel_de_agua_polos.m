%%TP1 Ejercicio 1, iniciso 2
pkg load control signal;
clc; clear all; close all;

%componentes
A1=1;
A2=1;
R1=1/2;
R2=1/3;
h1(1)=0;
h2(1)=0;

%matrices del espacio de estados (ss)
A=[-1/(A1*R1) 1/(A1*R1); 1/(A2*R1) -(1/(A2*R1)+1/(A2*R2))]
B=[1/A1; 0]
Ct=[0 1]
D=[0]
X=[h1 ; h2];
u(1)=0;
[num,den]=ss2tf(A,B,Ct,D);
sys=tf(num,den);

A_amp=[A zeros(2,1); -Ct 0]
B_amp=[B(:,1);0]

%%Diseño de la ubicacion de los polos
P= [-1.5+j -1.5-j -1]; %Polos deseados
K=place(A_amp,B_amp,P) %%Ganancia de realimentacion

Acl=A_amp-B_amp*K;
eig(Acl) %%Comprobamos que los plos estan donde queremos

%Kr=1/dcgain(ss(Acl,B_amp,[0 1 0],D)); %%Ganancia de entrada para corregir el offset


%Variables de utilidad
%%Para determinar h se hace lo siguiente:
f=max(abs(pole(sys)))%%Esto nos imprime la frecuencia del polo más alta del sistema.

%Como por el teorema del muestreo, tenemos que muestrear al doble de esa frecuencia
%determinamos h como 1/10f, donde f es la frecuencia más alta
h=1/(2*f)
tiempo=(10/h);
t=0:h:(tiempo*h);
i=1;
referencia=1; %%altura en metros del tanque 2
psi=0; %Error integrado
while(i<(tiempo+2))

  h1(i)=X(1);h2(i)=X(2);

  u=-K*([X;psi]);

  accion(i)=u;

  X_P=A*X+B*u;%X punto
  X=X+h*X_P;%Esto es el cálculo de la integral como sumatoria
  psi=psi+h*(referencia-Ct*X);

  i=i+1;

end

%Imprimo como varían mis variables de estado y mi entrada
figure 1;
subplot(3,1,1) ;plot(t,h1); title("altura del primer tanque"); grid on;
subplot(3,1,2) ;plot(t,h2); title("altura del segundo tanque");grid on; hold;
step(sys,"color","r"), legend("altura controlada","altura sin controlador");
subplot(3,1,3) ;plot(t,accion,"linewidth",2); title("caudal");grid on;

