## TP8: Control vía HBJ con estimación de estados
## Autor: Pedro Valentin, Nieva Miguel
## Consignas
## Implementar un controlador por la ecuación de HBJ con el estimador de Kalman
## para el caso del Avión

clear all; pkg load control;

##  Variables de simulacion
TamanioFuente=14;
randn('state',1);
Realizaciones=3; %Cantidad de realizaciones para el Monte Carlo.
Ts=0.01;

##----- Sigma determina la varianza del ruido -----##
sigma= 1e-5; colorc='k' %%Ruido nulo para no interferir con el sistema
##sigma= 0.05; colorc='r'
kmax=6000;
sQ=sigma; %Covarianza del ruido de estado Sigma=sqrt(sQ)
sR=sigma; %Covarianza del ruido de medicion sigma=sqrt(sR)
F_=sQ*eye(4);
G_=sR*eye(2);
I=eye(4);

##Sistema
a=0.07; w=9; b=5; c=150;

altura(1)=500;
x=[0;0;0;altura(1)];
x0=x;

alfa=zeros(Realizaciones,kmax);
fi=zeros(Realizaciones,kmax);
fi_p=zeros(Realizaciones,kmax);
altura=zeros(Realizaciones,kmax);

alfa(1)=x(1); fi(1)=x(2); fi_p(1)=x(3); altura(1)=x(4);

Mat_Ac=[-a a 0 0;
         0 0 1 0;
        w^2 -w^2 0 0;
         c 0 0 0]
Mat_Bc=[ 0;
         0;
         b*w^2;
         0]
Mat_C=[ 0 0 0 1;
        0 1 0 0];

sys_c=ss(Mat_Ac,Mat_Bc,Mat_C,0);
sys_d=c2d(sys_c,Ts,'zoh');

Mat_A=sys_d.a;
Mat_B=sys_d.b;
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B ];%Matriz Controlabilidad
rango=rank(Mat_M);

if rango==rank(Mat_A)
  disp("El sistema es controlable");
else
  disp("El sistema no es controlable");
  disp("El rango del sistema es");
  rank(Mat_A)
  disp("Pero la matriz de controlabilidad es");
  rango
end

%_____________ESTIMADOR KALMAN______________
P_Kalman=F_*F_'; #Po, ec 6-33
P11=zeros(1,5000);P22=P11;P33=P11;P44=P11;
%%La ganancia de Kalman es representada por "T" en el apunte
for h_k=1:5000
  P_Kalman_=Mat_A*P_Kalman*Mat_A'+(F_*F_'); #Ec 6-38, es Pk⁻
  K_Kalman=(P_Kalman_*Mat_C')/(Mat_C*P_Kalman_*Mat_C'+(G_*G_')); %Ganancia de Kalman, ec 6-31
  P_Kalman=(eye(4)-K_Kalman*Mat_C)*P_Kalman_; #Ec 6-35

  ##Asignacion de las Covarianzas de la matriz de covarianzas
  P11(h_k)=P_Kalman(1,1);
  P22(h_k)=P_Kalman(2,2);
  P33(h_k)=P_Kalman(3,3);
  P44(h_k)=P_Kalman(4,4);
end
EK=abs(eig(Mat_A-K_Kalman*Mat_C));


%_____________ Método HJB ______________
c_ai= poly(eig(Mat_Ac));

Mat_W=[c_ai(4) c_ai(3) c_ai(2) 1;c_ai(3) c_ai(2) 1 0;c_ai(2) 1 0 0;1 0 0 0];
Mat_M=[Mat_Bc Mat_Ac*Mat_Bc Mat_Ac^2*Mat_Bc Mat_Ac^3*Mat_Bc ];%Matriz Controlabilidad continua

Mat_T=Mat_M*Mat_W;

A_controlable=inv(Mat_T)*Mat_Ac*Mat_T; %Verificación de que T esté bien

a4=-A_controlable(4,1);%a4
a3=-A_controlable(4,2);
a2=-A_controlable(4,3);
a1=-A_controlable(4,4);

q1=1e1; #Dirección de Vuelo
q2=1e2; #Ángulo de Inclinacion
q3=2e-2; #Velocidad de Inclinacion
q4=1e1; #Altura
R=1e0;

p1=.5*(-4*a4*R+sqrt((4*a4*R)^2+16*q1*R)); %%Ec 7-80
p2=.5*(-4*a3*R+sqrt((4*a3*R)^2+16*q2*R)); %%Ec 7-81
p3=.5*(-4*a2*R+sqrt((4*a2*R)^2+16*q3*R)); %%Ec 7-82
p4=.5*(-3*a1*R+sqrt((3*a1*R)^2+8*q4*R)); %%Ec 7-85

K=(([p1 p2 p3 p4])/(2*R))*inv(Mat_T); %%Ec 7-88


if any(real(eig(Mat_Ac-Mat_Bc*K))>0) %%Autovalores del sistema realimentado
  eig(Mat_Ac-Mat_Bc*K)
  error("SISTEMA INESTABLE");
else
  eig(Mat_Ac-Mat_Bc*K)
endif
##----- Simulación -----##
u=zeros(Realizaciones,kmax);
y_sal=zeros(Realizaciones,kmax);
t=0:kmax-1;

for trial=1:Realizaciones %Empieza el Monte Carlo
  v=randn(4,kmax);%Señales aleatorios de media nula y varianza unidad.
  w=randn(2,kmax);
  x=x0+F_*v(:,1);
  alfa(trial,1)=x0(1);
  fi(trial,1)=x0(2);
  fi_p(trial,1)=x0(3);
  altura(trial,1)=x0(4);
  x_hat=zeros(size(x0));
  x_hat_=x_hat;

  for ki=1:kmax-1
    Y=Mat_C*x+G_*w(:,ki);

    x_hat=x_hat_+K_Kalman*(Y-Mat_C*x_hat_); #Estimador de Kalman

    u(trial,ki)=-K*x_hat; #HJB del sistema
    u(trial,ki) = max(min(u(trial,ki), 0.3), -0.3); #Se limita la acción de control

    x=modavion(Ts,x,u(trial,ki))+F_*v(:,ki);

    x_hat_=Mat_A*x_hat+Mat_B*u(trial,ki); #Estimador de Kalman


    alfa(trial,ki+1)=x(1);
    fi(trial,ki+1)=x(2);
    fi_p(trial,ki+1)=x(3);
    altura(trial,ki+1)=x(4);

  end

end

t=t*Ts;


##----- Ploteos -----##
figure(1);hold on;
subplot(3,1,1);
hold on;grid on; title('Direccion de vuelo','FontSize',TamanioFuente-1);
plot(t,mean(alfa),colorc);
plot(t,mean(alfa)+.5*sqrt(var(alfa)),["--",colorc]);
plot(t,mean(alfa)- .5*sqrt(var(alfa)),["--",colorc]);
##
##subplot(3,2,2);
##hold on;grid on;title('Inclinacion del avion','FontSize',TamanioFuente-1);
##plot(t,mean(fi),colorc);
##plot(t,mean(fi)+.5*sqrt(var(fi)),["--",colorc]);
##plot(t,mean(fi)-.5*sqrt(var(fi)),["--",colorc]);

##subplot(3,2,3);
##grid on;title('Velocidad de inclinacion','FontSize',TamanioFuente-1);hold on;
##plot(t,mean(fi_p),colorc);
##plot(t,mean(fi_p)+.5*sqrt(var(fi_p)),["--",colorc]);
##plot(t,mean(fi_p)- .5*sqrt(var(fi_p)),["--",colorc]);

subplot(3,1,2);
hold on; grid on;title('Altura','FontSize',TamanioFuente-1);hold on;
plot(t,mean(altura),colorc);
plot(t,mean(altura)+.5*sqrt(var(altura)),["--",colorc]);
plot(t,mean(altura)- .5*sqrt(var(altura)),["--",colorc]);

subplot(3,1,3);
grid on;title('Acción de control','FontSize',TamanioFuente-1);xlabel('Tiempo en Seg.','FontSize',TamanioFuente-1);hold on;
plot(t,mean(u),colorc);
plot(t,mean(u)+.5*sqrt(var(u)),["--",colorc]);
plot(t,mean(u)-.5*sqrt(var(u)),["--",colorc]);

