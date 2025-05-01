pkg load control signal;
# Tiene formato de código
%% Ajuste por mínimos cuadrados. Algoritmo que se emplea luego del analisis
% de correlación y densidad espectral de potencia.
% Éste algoritmo funciona porque se ha ajustado a la señal de entrada y al
% tiempo de muestro a valores adecuados.
clear all; close all;clc;
Med=2500;
ts=1/50;
orden_b=2; %Orden del Numerador
orden_a=3; %Orden del Denominador
t=0:ts:20000*ts;
randn('state',0);

% Señal de entrada
StepAmplitude=1;
ue=StepAmplitude*sign(sin(2*pi*.2*1.0*t));
m=7;
x=ones(m,1);
N=length(ue);%Puntos de la PRBS para muestrear
el=50;%Repeticiones del ZOH
for k=1:el:N
    n_b=xor(x(7),x(6));
    y(k:k+el-1)=x(7);
    for h=m-1:-1:1
        x(h+1)=x(h);
    end
    x(1)=n_b; %Ingreso el nuevo valor
end
% x=[]; % plot(y) % break
ue=y(1:N); %Uso de PRBS

%Planta
num=[2 -1];
den=[1 3 6];
sys=tf(num,den);
[y_D,t_D]=lsim(sys,ue,t);
ys=y_D';

%Cuanto es lo que tarda en iniciar el algoritmo
off_set=orden_a+1+50;
u=zeros(Med,1);
z=zeros(Med,1);
u=ue(1+off_set:length(u)+off_set)';ui=u;
z=ys(1+off_set:length(z)+off_set)';zi=z;
for jj=orden_a+1:Med
    vec_a=fliplr([ u(jj-orden_b:jj-1); -z(jj-orden_a:jj-1)]');
    H(jj-orden_a,:)=vec_a;
end
[aa bb ]=size(H);
Z=(z(orden_a+1:end));
in_1=H';
in_2=in_1*H;
in_3=inv(in_2);
in_4=in_3*in_1;
c=in_4*Z; % Ec 131
abs(roots([1; c(1:orden_a)]));
zo=z;
z=zeros(Med,1);u=ui;
z=zi(1:orden_a);
u=ue(1+off_set:off_set+Med)';
for k=orden_a+1:length(u)
    zt=-flip(z(k-orden_a:k-1));
    ut=flip(u(k-orden_b:k-1));
    z(k)=c'*[zt;ut];    %Ec 125
end
dend=[1; c(1:orden_a)]';numd=[c(orden_a+1:end)]';
sys_id=tf(numd,dend,ts,'Variable','z^-1');
%ue=sign(sin(2*pi*.010*t));

[y_sal,t_sal]=lsim(sys_id,ue,t);
[y_D,t_D]=lsim(sys,ue,t);

subplot(2,1,1);hold on;
plot(z,'.'); grid on;
plot(zo,'+r');
plot(ys(1+off_set:off_set+Med),'y');
legend('Estimada','Mediciones','Real');
title(['Ajuste con el orden b_s/a_s:' num2str(orden_b) '/' num2str(orden_a)]);
xlabel('Muestras')

subplot(2,1,2)
hold on; grid on;
plot(t_D*ts,y_D,'.');
plot(t_sal*ts,y_sal,'k');legend('Real','Identificada')
title('Desempeño del modelo ajustado');xlabel('Tiempo. [Seg.]')

%%Respuesta al escalon lo calculamos así para luego exportarla
t2=0:0.02:5;
[y_D2,t_D2]=lsim(sys,ones(1,length(t2)).*StepAmplitude,t2);
[y_sal2,t_sal2]=lsim(sys_id,ones(1,length(t2)).*StepAmplitude,t2);

figure
hold on; grid on;
plot(t_D2,y_D2,'.');
plot(t_sal2,y_sal2,'k');
legend('Real','Identificada');
title('Step Response');xlabel('Tiempo. [Seg.]')

% sys %Sistema original
sys_Norm=sys;
sys_Norm.num{1}=sys.num{1}/sys.den{1}(end);
sys_Norm.den{1}=sys.den{1}/sys.den{1}(end);

sysc = d2c(sys_id,'tustin'); %sistema identificado
sysc_Norm=tf(sysc.num{1}/sysc.den{1}(end),sysc.den{1}/sysc.den{1}(end));

disp("El sistema continuo original normalizado es:")
sys_Norm
disp("El sistema continuo identificado es:")
sysc_Norm

%%Realizamos el bode de esta forma para exportar luego los resultados a python
[mag_sys,fase_sys,W]=bode(sys,logspace(-1,1));
[mag_sysid,fase_sysid,Wid]=bode(sys_id,logspace(-1,1));


figure
subplot(2,1,1); grid on; hold on;
semilogx(W,mag_sys); %Bode real
semilogx(Wid,mag_sysid); %Bode Identificado
legend('Original','Identificada');

subplot(2,1,2); grid on; hold on;
semilogx(W,fase_sys);
semilogx(Wid,fase_sysid);

save -v7 TP5_caso_a.mat
