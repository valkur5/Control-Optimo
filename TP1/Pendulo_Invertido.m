%%TP1 Ejercicio 2, iniciso 1
pkg load control signal symbolic;
clc; clear all; close all;

%Parámetros
m  = 0.1;
F  = 0.1;
l  = 1.6;
g  = 9.8;
M  = 1.5;

Ts=0.01;
h = 1e-4;
t = 0:h:55;
p_max = floor(7/Ts);
%Matrices de estado
A=[0       1               0        0;
      0    -F/M            -m*g/M      0;
      0       0               0        1;
      0   F/(l*M)     g*(m+M)/(l*M)  0];

B=[   0   ;
      1/M  ;
       0   ;
   -1/(l*M)];

C = [1 0 0 0];

D=0;

%Sistema discreto
sysDisc = c2d(ss(A,B,C,D), Ts, 'zoh');
Ad     = sysDisc.a;
Bd     = sysDisc.b;
Cd     = sysDisc.c;
Dd     = sysDisc.d;

% Sistema ampliado
Aa = [Ad , zeros(4,1) ; -Cd(1,:)*Ad, eye(1)];
Ba = [Bd; -Cd(1,:)*Bd];

Qamp = diag([100 100 100 10 1]);
Ramp = diag([1]);

K_c = dlqr(Aa,Ba,Qamp,Ramp);
%%Calculamos ahora el observador
Ao = Ad';
Bo = Cd';
Co = Bd';

Qo = diag([10 10 10 10]);
Ro = diag([1]);
Ko = (dlqr(Ao,Bo,Qo,Ro))';

%%Condiciones de simulación
d(1)     = 0;
d_p(1)   = 0;
phi(1)   = 0.2; %Posición inicial
phi_p(1) = 0;
phi_pp(1) = 0;

ref = 1; %Posición final

X= [d(1) d_p(1) phi(1) phi_p(1)]';


u(1) = 0;
ei(1) = 0;
x_hat = [0 0 0 0]';

ii=1;
%Simulación
for ki=1:p_max

    % Salida de dos componentes
    Ys   = Cd*X;             % salida del sistema
    Y_obs = Cd*(x_hat);   % salida del observador

    ei(ki+1)= ei(ki)+ref-Ys(1);

    %Ley de control
    u1(ki)  = -K_c*[(x_hat);ei(ki+1)];     % con observador


    % Integraciones de Euler por paso de simulaci�n
    for kii=1:Ts/h
        u(ii) = u1(ki);

        % C�lculo por sistema no lineal
        d_pp       = (1/(M+m))*(u(ii) - m*l*phi_pp*cos(phi(ii)) + m*l*phi_p(ii)^2*sin(phi(ii)) - F*d_p(ii));
        phi_pp     = (1/l)*(g*sin(phi(ii)) - d_pp*cos(phi(ii)));
        d_p(ii+1)   = d_p(ii) + h*d_pp;
        d(ii+1)     = d(ii) + h*d_p(ii);
        phi_p(ii+1) = phi_p(ii) + h*phi_pp;
        phi(ii+1)   = phi(ii) + h*phi_p(ii);

        ii=ii+1;
    end

    % Actualizaci�n de los estados
    % Estados del sistema
    X     = [d(ii-1) d_p(ii-1) phi(ii-1) phi_p(ii-1)]';
    % Estados estimados por el observador
    x_hat = Ad*(x_hat) + Bd*u1(ki) + Ko*(Ys - Y_obs);
end
u(ii) = u1(ki);


subplot(2,1,1); plot(phi);grid on; title("phi")
subplot(2,1,2);plot(d);grid on; title("d")
