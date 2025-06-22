## TP8: Control vía HBJ con estimación de estados
## Autor: Pedro Valentin, Nieva Miguel
## Consignas
## Implementar un controlador por la ecuación de HBJ con el estimador de Kalman
## para el caso del péndulo invertido
## Luego intentar estabilizar en equilibrio inestable partiendo de +-pi
clear all; pkg load control;

##  Variables de simulacion
TamanioFuente=14;
randn('state',1);
Realizaciones=1; %Cantidad de realizaciones para el Monte Carlo.
Ts=0.01;

##----- Por consigna, Sigma debe ser 0, 0.01, 0.02, 0.05, 0.1 -----##
sigma= 0.01
kmax=2000;
sQ=sigma; %Covarianza del ruido de estado Sigma=sqrt(sQ)
sR=sigma; %Covarianza del ruido de medicion sigma=sqrt(sR)
F_=sQ*eye(4);
G_=sR*eye(2);
I=eye(4);

##Sistema
m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;

alfa(1)=0.8;
x=[0;0;alfa(1);0];
x0=x;
p(1)=x(1); p_p(1)=x(2); alfa(1)=x(3); omega(1)=x(4);
Mat_Ac=[0 1 0 0;
        0 -Fricc/M -m*g/M 0;
        0 0 0 1;
        0 Fricc/(long*M) g*(m+M)/(long*M) 0];
Mat_Bc=[0;
        1/M;
        0;
        -1/(long*M)];
Mat_C=[1 0 0 0;
       0 0 1 0];

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

%_____________Método HJB ______________
c_ai= poly(eig(Mat_Ac));

Mat_W=[c_ai(4) c_ai(3) c_ai(2) 1;c_ai(3) c_ai(2) 1 0;c_ai(2) 1 0 0;1 0 0 0];
Mat_M=[Mat_Bc Mat_Ac*Mat_Bc Mat_Ac^2*Mat_Bc Mat_Ac^3*Mat_Bc ];%Matriz Controlabilidad continua

Mat_T=Mat_M*Mat_W;

A_controlable=inv(Mat_T)*Mat_Ac*Mat_T; %Verificación de que T esté bien

a4=-A_controlable(4,1);%a4
a3=-A_controlable(4,2);
a2=-A_controlable(4,3);
a1=-A_controlable(4,4);

q1=1e1;q2=1e0;q3=1e4;q4=1e4;
R=1e0;

p1=.5*(-4*a4*R+sqrt((4*a4*R)^2+16*q1*R)); %%Ec 7-80
p2=.5*(-4*a3*R+sqrt((4*a3*R)^2+16*q2*R)); %%Ec 7-81
p3=.5*(-4*a2*R+sqrt((4*a2*R)^2+16*q3*R)); %%Ec 7-82
p4=.5*(-3*a1*R+sqrt((3*a1*R)^2+8*q4*R)); %%Ec 7-85

K=(([p1 p2 p3 p4])/(2*R))*inv(Mat_T); %%Ec 7-88

eig(Mat_Ac-Mat_Bc*K) %%Autovalores del sistema realimentado

##----- Simulación -----##
u=zeros(Realizaciones,kmax);
y_sal=zeros(Realizaciones,kmax);
t=0:kmax-1;

for trial=1:Realizaciones %Empieza el Monte Carlo
  v=randn(4,kmax);%Señales aleatorios de media nula y varianza unidad.
  w=randn(2,kmax);
  x=x0+F_*v(:,1);
  p(trial,1)=x0(1);
  p_p(trial,1)=x0(2);
  alfa(trial,1)=x0(3);
  omega(trial,1)=x0(4);
  x_hat=zeros(size(x0));
  x_hat_=x_hat;

  for ki=1:kmax-1
    Y=Mat_C*x+G_*w(:,ki);
    y_sal(trial,ki)=Y(1);
    x_hat=x_hat_+K_Kalman*(Y -Mat_C*x_hat_); #Estimador de Kalman
    u(trial,ki)=-K*x_hat;#LQG del sistema
    x=mopdm(Ts,x,u(trial,ki))+F_*v(:,ki);

    Y_O=Mat_C*x_hat;
    y_sal_O(trial,ki)=Y_O(1);

    x_hat_=Mat_A*x_hat+Mat_B*u(trial,ki); colorc='r'; Ley="Kalman"; #Estimador de Kalman


    p(trial,ki+1)=x(1);
    p_p(trial,ki+1)=x(2);
    alfa(trial,ki+1)=x(3);
    omega(trial,ki+1)=x(4);

  end

end

t=t*Ts;


##----- Ploteos -----##
figure(1);hold on;
subplot(3,2,1);
hold on;grid on; title('Velocidad ángulo','FontSize',TamanioFuente-1);
plot(t,mean(omega),colorc);
plot(t,mean(omega)+.5*sqrt(var(omega)),["--",colorc]);
plot(t,mean(omega)-.5*sqrt(var(omega)),["--",colorc]);
##xlim([0 0.45])

subplot(3,2,2);
hold on;grid on;title('Ángulo','FontSize',TamanioFuente-1);
plot(t,mean(alfa),colorc);
plot(t,mean(alfa)+.5*sqrt(var(alfa)),["--",colorc]);
plot(t,mean(alfa)-.5*sqrt(var(alfa)),["--",colorc]);
##xlim([0 2])

subplot(3,2,3);
grid on;title('Posición carro','FontSize',TamanioFuente-1);hold on;
plot(t,mean(p),colorc);
plot(t,mean(p)+.5*sqrt(var(p)),["--",colorc]);
plot(t,mean(p)-.5*sqrt(var(p)),["--",colorc]);

subplot(3,2,4);
hold on; grid on;title('Velocidad carro','FontSize',TamanioFuente-1);hold on;
plot(t,mean(p_p),colorc);
plot(t,mean(p_p)+.5*sqrt(var(p_p)),["--",colorc]);
plot(t,mean(p_p)-.5*sqrt(var(p_p)),["--",colorc]);
##xlim([0 2])

subplot(3,1,3);
grid on;title('Acción de control','FontSize',TamanioFuente-1);xlabel('Tiempo en Seg.','FontSize',TamanioFuente-1);hold on;
plot(t,mean(u),colorc);
plot(t,mean(u)+.5*sqrt(var(u)),["--",colorc]);
plot(t,mean(u)-.5*sqrt(var(u)),["--",colorc]);
##xlim([0 0.45])
##axes( 'visible', 'off', 'title', ['Pendulo Invertido, Angulo inicial \Phi= ',num2str(alfa(1))],"FontSize",TamanioFuente );
