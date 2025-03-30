%%TP1 Ejercicio 1, iniciso 1
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
sys=ss(A,B,Ct,D);
Y(1)=0;

%Variables de utilidad
%%Para determinar h se hace lo siguiente:
f=max(abs(pole(sys)))%%Esto nos imprime la frecuencia del polo más alta del sistema.

%Como por el teorema del muestreo, tenemos que muestrear al doble de esa frecuencia
%determinamos h como 1/10f, donde f es la frecuencia más alta
h=1/(2*f) %El h como yo
tiempo=(10/h);
t=0:h:(tiempo*h);
i=1; k=0;
referencia=1; %%altura en metros del tanque 2
%%PID
Kp=3; Ki=3;Kd=0.01;

A1=((2*Kp*h)+(Ki*(h^2))+(2*Kd))/(2*h);
B1=(-2*Kp*h+Ki*(h^2)-4*Kd)/(2*h);
C1=Kd/h;
error(1)=0;


while(i<(tiempo+2))
  k=i+2;
  h1(i)=X(1);h2(i)=X(2);
  error(k)=referencia-h2(i); %%calculamos el error

  u=[u+A1*error(k)+B1*error(k-1)+C1*error(k-2)]; %%calculamos la acción de control

  accion(i)=u;

  for ii=1:10
      X_P=A*X+B*u;%X punto
      X=X+.1*h*X_P;%Esto es el cálculo de la integral como sumatoria

  endfor


  i=i+1;

end

%Imprimo como varían mis variables de estado y mi entrada
figure 1;
subplot(3,1,1) ;plot(t,h1); title("altura del primer tanque"); grid on;
subplot(3,1,2) ;plot(t,h2); grid on; hold; step(sys,t,"color","r"), legend("altura controlada","altura sin controlador"); title("altura del segundo tanque");
subplot(3,1,3) ;plot(t,accion); title("caudal u(t)");grid on;

