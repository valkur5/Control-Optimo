%%TP1 Ejercicio 1, iniciso 4
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
Ts=0.04;
h=Ts/20;
muestras=1000;
t = 0:h:(muestras*Ts);
i=1; k=0;
referencia=1; %%altura en metros del tanque 2


%%DLQR
%%Primero discretizamos el sistema como:
sys_D=c2d(sys,h,'zoh');
Ad  = sys_D.a;
Bd  = sys_D.b;
Ctd = sys_D.c;
Dd  = sys_D.d;

A_ampd=[Ad zeros(2,1); -Ctd*Ad 1]
B_ampd=[Bd;-Ctd*Bd]

Qd  = diag([10 .1 100])
Rd   = 1
Kd=dlqr(A_ampd,B_ampd,Qd,Rd);
psi=0;
Y(1)=0;
u1(1)=0;

for ki=1:muestras
  u1(ki)=-Kd*([X;psi]);
  for kii=1:Ts/h
    u(i)=u1(ki);

    h2_p=(-(h2(i)-h1(i))/A2*R1)-(h2(i)/R2*A2);
    h1_p=(h2(i)-h1(i))/A1*R1+u(i)/A1;
    h2(i+1)=h2(i)+h*h2_p;
    h1(i+1)=h1(i)+h*h1_p;
    psi=psi+h*(referencia-h2(i));
    X(1)=h1(i);X(2)=h2(i);
    i=i+1;
  endfor

endfor
u(i)=u1(ki);
%Imprimo como varían mis variables de estado y mi entrada
figure 1;
subplot(3,1,1) ;plot(t,h1); title("altura del primer tanque"); grid on;
subplot(3,1,2) ;plot(t,h2); grid on; hold;
step(sys,t,"color","r"), legend("altura controlada","altura sin controlador"); title("altura del segundo tanque");
subplot(3,1,3) ;plot(t,u); title("caudal");grid on;

