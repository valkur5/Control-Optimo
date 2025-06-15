clear all; pkg load control;
randn('state',1);

TamanioFuente=14;
%Condiciones iniciales
altura(1)=500; color='.r';colorc='r';
% alfa(1)=.5; color='.g';colorc='g';
% alfa(1)=.8; color='.b';colorc='b';
Realizaciones=50; %Cantidad de realizaciones para el Monte Carlo.

%Versión linealizada en el equilibrio inestable. Sontag Pp 104.
a=0.07;
w=9;
b=5;
c=150;

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

##----- Por consigna, Sigma debe ser 0, 0.01, 0.02, 0.05, 0.1
sigma= 0.01
kmax=5000;
sQ=sigma; %Covarianza del ruido de estado Sigma=sqrt(sQ)
sR=sigma; %Covarianza del ruido de medicion sigma=sqrt(sR)

F_=sQ*eye(4);
G_=sR*eye(2);
I=eye(4);

sys_c=ss(Mat_Ac,Mat_Bc,Mat_C,0);
Ts=0.01;
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

x=[0;0;0;altura(1)];
x0=x;
alfa(1)=x(1); fi(1)=x(2); fi_p(1)=x(3); altura(1)=x(4);
Aa=Mat_A;Ba=Mat_B;

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

##figure(4);hold on;
##semilogy(P11,"b");
##semilogy(P22,'g');
##semilogy(P33,'c');
##semilogy(P44,'r');
##legend("P11","P22","P33","P44");
##title('Evolución de P_1_1,P_2_2,P_3_3y P_4_4.')
##xlabel('Iteraciones');

EK=abs(eig(Mat_A-K_Kalman*Mat_C));
Q=diag([5 5 5 1e-4]);%Matrices de diseño del controlador DLQG
S=Q;
P=S; %condición inicial de P
R=1e0;

##----- DLQG para Kalman -----##
Kk=zeros(kmax,4);
for hi=kmax-1:-1:1
  P= Q + Aa'*P*Aa - Aa'*P*Ba/(R+Ba'*P*Ba) * Ba'*P*Aa; ##Ec 6-44 de la clase 6
  Kk(hi,:)=(R+Ba'*P*Ba)\Ba'*P; ## Ec 6-47 de la clase 6
  Ea(:,hi)=eig(Aa-Ba*Kk(hi,:)*Aa);
end

#plot(abs(Ea)')

##----- DLQR para la planta -----##
[Klqr,Pa,Elqr]=dlqr(Aa,Ba,Q,R);

##----- Cálculo del observador de Luenberger -----##
Qo=diag([1000 2000 10000 1000000000]);
Ro=diag([1e5, 1e-4]);
Aad=Mat_A';
Bad=Mat_C';

[Ko,Po,Em]=dlqr(Aad,Bad,Qo,Ro);
Ko=Ko';


Jmin=x0'*P*x0;J=0;t=0:kmax-1;u=zeros(Realizaciones,kmax);Jn_=zeros(Realizaciones,kmax);
y_sal=zeros(Realizaciones,kmax);
tic
for trial=1:Realizaciones %Empieza el Monte Carlo
  v=randn(4,kmax);%Señales aleatorios de media nula y varianza unidad.
  w=randn(2,kmax);
  x=x0+F_*v(:,1);
  alfa(trial,1)=x(1);
  fi(trial,1)=x(2);
  fi_p(trial,1)=x(3);
  altura(trial,1)=x(4);
  x_hat=zeros(size(x0));
  x_hat_=x_hat;
  for ki=1:kmax-1
    Y=Mat_C*x+G_*w(:,ki+1);
    y_sal(trial,ki)=Y(1);

    x_hat=x_hat_+K_Kalman*(Y-Mat_C*x_hat_); #Estimador de Kalman
    %u(trial,ki)=-Klqr*x_hat;
    u(trial,ki)=-Kk(ki,:)*x_hat; #DLQG del sistema

    x=modavion(Ts,x,u(trial,ki))+F_*v(:,ki+1);
    Y_O=Mat_C*x_hat;
    y_sal_O(trial,ki)=Y_O(1);
    Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*eye(4)*x + u(trial,ki)'*1*u(trial,ki));
    x_hat_=Mat_A*x_hat+Mat_B*u(trial,ki);colorc='k'; Ley="Kalman"; #Estimador de Kalman

##    x_hat=Mat_A*x_hat+Mat_B*u(trial,ki)+Ko*(Y-Mat_C*x_hat);colorc="r"; Ley="Luenberger"; #Observador de Luenberger

    alfa(trial,ki+1)=x(1);
    fi(trial,ki+1)=x(2);
    fi_p(trial,ki+1)=x(3);
    altura(trial,ki+1)=x(4);

##    P_Kalman_=Mat_A*P_Kalman*Mat_A'+(F_*F_'); #Ec 6-38, es Pk⁻
##    K_Kalman=(P_Kalman_*Mat_C')/(Mat_C*P_Kalman_*Mat_C'+(G_*G_')); %Ganancia de Kalman, ec 6-31
##    P_Kalman=(eye(4)-K_Kalman*Mat_C)*P_Kalman_; #Ec 6-35
  end
  Jn_(trial,ki+1)=Jn_(trial,ki+1)+x'*S*x;
end
toc
t=t*Ts;Jn=mean(Jn_);
disp(["Utilizando " Ley ' Jn(end)=' num2str(Jn(end)) '. Altura inicial = ' num2str(altura(1)) '[mts].']);

figure(1);hold on;
subplot(3,2,1);hold on;grid on; title('Direccion de vuelo','FontSize',TamanioFuente-1);hold on;
plot(t,mean(alfa),colorc); hold on;
plot(t,mean(alfa)+.5*sqrt(var(alfa)),["--",colorc]);
plot(t,mean(alfa)- .5*sqrt(var(alfa)),["--",colorc]);

subplot(3,2,2);hold on;grid on;title('Inclinacion del avion','FontSize',TamanioFuente-1);hold on;
plot(t,mean(fi),colorc); hold on;
plot(t,mean(fi)+.5*sqrt(var(fi)),["--",colorc]);
plot(t,mean(fi)-.5*sqrt(var(fi)),["--",colorc]);

subplot(3,2,3);hold on; grid on;title('Velocidad de inclinacion','FontSize',TamanioFuente-1);hold on;
plot(t,mean(fi_p),colorc); hold on;
plot(t,mean(fi_p)+.5*sqrt(var(fi_p)),["--",colorc]);
plot(t,mean(fi_p)- .5*sqrt(var(fi_p)),["--",colorc]);

subplot(3,2,4);hold on; grid on;title('Altura','FontSize',TamanioFuente-1);hold on;
plot(t,mean(altura),colorc); hold on;
plot(t,mean(altura)+.5*sqrt(var(altura)),["--",colorc]);
plot(t,mean(altura)- .5*sqrt(var(altura)),["--",colorc]);

subplot(3,1,3); grid on;title('Acción de control','FontSize',TamanioFuente-1);
xlabel('Tiempo en Seg.','FontSize',TamanioFuente-1);hold on;
plot(t,mean(u),colorc); hold on;
plot(t,mean(u)+.5*sqrt(var(u)),["--",colorc]);
plot(t,mean(u)- .5*sqrt(var(u)),["--",colorc]);

axes( 'visible', 'off', 'title', ['Avion, Altura inicial h= ',num2str(altura(1))],"FontSize",TamanioFuente );
#save -v7 TP7_Avion_Luenberger.mat
