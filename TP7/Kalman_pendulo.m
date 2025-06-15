clear all; pkg load control;
randn('state',1);

m=.1;Fricc=0.1; long=0.6;g=9.8;M=.5;
TamanioFuente=14;
%Condiciones iniciales
alfa(1)=.8;
% alfa(1)=.5; color='.g';colorc='g';
% alfa(1)=.8; color='.b';colorc='b';
Realizaciones=50; %Cantidad de realizaciones para el Monte Carlo.

%Versión linealizada en el equilibrio inestable. Sontag Pp 104.

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

##----- Por consigna, Sigma debe ser 0, 0.01, 0.02, 0.05, 0.1 -----##
sigma= 0.1
kmax=2000;
sQ=sigma; %Covarianza del ruido de estado Sigma=sqrt(sQ)
sR=sigma; %Covarianza del ruido de medicion sigma=sqrt(sR)

F_=sQ*eye(4);
G_=sR*eye(2);
I=eye(4);

sys_c=ss(Mat_Ac,Mat_Bc,Mat_C,0);Ts=0.01;
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

x=[0;0;alfa(1);0];x0=x;
p(1)=x(1); p_p(1)=x(2); alfa(1)=x(3); omega(1)=x(4);
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
Q=1e3*diag([1e1 1e1 1e1 1e1]);%Matrices de diseño del controlador DLQG
Q=diag([1e2 1e1 1e6 1e1]); #Matriz obtenida en el TP6 del pendulo
S=Q;
P=S; %condición inicial de P
R=1e5;
R=1e0; #Matriz obtenida en el TP6
Kk=zeros(kmax,4);
for hi=kmax-1:-1:1
  P= Q + Aa'*P*Aa - Aa'*P*Ba/(R+Ba'*P*Ba) * Ba'*P*Aa; ##Ec 6-44 de la clase 6
  Kk(hi,:)=(R+Ba'*P*Ba)\Ba'*P; ## Ec 6-47 de la clase 6
  Kv(hi,:)=inv(R+Ba'*P*Ba)*Ba'*P*F_;
  Ea(:,hi)=eig(Aa-Ba*Kk(hi,:)*Aa);
end

#plot(abs(Ea)')

[Klqr,Pa,Elqr]=dlqr(Aa,Ba,Q,R);

%Cálculo del observador
Qo=diag([1e0 1e0 1e0 1e0]);
Ro=diag([1e4 1e0]);
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
##    u(trial,ki)=-Klqr*x_hat; #LQR del sistema
    u(trial,ki)=-Kk(ki,:)*x_hat;#LQG del sistema
    x=mopdm(Ts,x,u(trial,ki))+F_*v(:,ki);
    Y_O=Mat_C*x_hat;
    y_sal_O(trial,ki)=Y_O(1);
    Jn_(trial,ki+1)=Jn_(trial,ki)+(x'*eye(4)*x + u(trial,ki)'*1*u(trial,ki));

    x_hat_=Mat_A*x_hat+Mat_B*u(trial,ki); colorc='r'; Ley="Kalman"; #Estimador de Kalman

##    x_hat=Mat_A*x_hat+Mat_B*u(trial,ki)+Ko*(Y-Mat_C*x_hat); colorc='r'; Ley="Luenberger"; #Observador de Luenberger

    p(trial,ki+1)=x(1);
    p_p(trial,ki+1)=x(2);
    alfa(trial,ki+1)=x(3);
    omega(trial,ki+1)=x(4);

    if ki==round((kmax-1)/3)
      F_=2*sQ*eye(4);
      G_=2*sR*eye(2);
##      for h_k=1:5000
##        P_Kalman_=Mat_A*P_Kalman*Mat_A'+(F_*F_'); #Ec 6-38, es Pk⁻
##        K_Kalman=(P_Kalman_*Mat_C')/(Mat_C*P_Kalman_*Mat_C'+(G_*G_')); %Ganancia de Kalman, ec 6-31
##        P_Kalman=(eye(4)-K_Kalman*Mat_C)*P_Kalman_; #Ec 6-35
##      endfor
    elseif ki==round((kmax-1)*2/3)
      F_=sQ*eye(4);
      G_=sR*eye(2);
##      for h_k=1:5000
##        P_Kalman_=Mat_A*P_Kalman*Mat_A'+(F_*F_'); #Ec 6-38, es Pk⁻
##        K_Kalman=(P_Kalman_*Mat_C')/(Mat_C*P_Kalman_*Mat_C'+(G_*G_')); %Ganancia de Kalman, ec 6-31
##        P_Kalman=(eye(4)-K_Kalman*Mat_C)*P_Kalman_; #Ec 6-35
##      endfor
    endif
##      P_Kalman_=Mat_A*P_Kalman*Mat_A'+(F_*F_'); #Ec 6-38, es Pk⁻
##      K_Kalman=(P_Kalman_*Mat_C')/(Mat_C*P_Kalman_*Mat_C'+(G_*G_')); %Ganancia de Kalman, ec 6-31
##      P_Kalman=(eye(4)-K_Kalman*Mat_C)*P_Kalman_; #Ec 6-35
  end

  Jn_(trial,ki+1)=Jn_(trial,ki+1)+x'*S*x;
end
toc
t=t*Ts;Jn=mean(Jn_);
disp(["Utilizando " Ley ' Jn(end)=' num2str(Jn(end)) '. Inclinacion inicial = ' num2str(alfa(1)) '[rads].']);

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
axes( 'visible', 'off', 'title', ['Pendulo Invertido, Angulo inicial \Phi= ',num2str(alfa(1))],"FontSize",TamanioFuente );

#save -v7 TP7_2_Kalman_pendulo_Estacionario.mat
