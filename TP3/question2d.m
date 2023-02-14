%Script: Question 2.d du TP3.
clc
format compact;
format short e;

% Conditions initiales
t0=0;

g = 9.81;
theta = ( 45 + 31 / 60 ) / 360 * 2 *pi;
l=8;
Omega = 7.29*10^(-5);
omega = sqrt(g/l);
omegaZero = sqrt(omega^2+Omega^2*sin(theta)^2);

x0=0.5;
x0Prime = 0;
y0 = 0;
y0Prime = 0;
z0 = x0 + 1i * y0;
z0Prime = x0Prime + 1i * y0Prime;
z0Total = [z0,z0Prime];

%Temps finale
tf = 20*pi/omegaZero;

% Solution de référence
zRef=solexFoucault(tf,z0,z0Prime,omega,Omega,theta)

% Nombre de pas
n=2.^(14);
taille = length(n);

% Initialisation de vecteur
zf = zeros(length(taille));

% Méthode prédicteur correcteur
for i=1:taille;
ni=n(i); % Nombre de pas nécessaire pour cet itération

% Trajectoire
[ti,zi]=predcor(@foucaultODE,[t0,tf],z0Total,ni,'AB3','AM4',omega,Omega,theta); 
%[ti,zi]=predcor2(@foucaultODE,[t0,tf],z0Total,ni,'eulerunpas2','cnunpas2',omega,Omega,theta); 

zf(i)=zi(end,1);
end;
zf
anglePosition = asin(imag(zf)/real(zf))
nombreOscillations = (2*pi/anglePosition)*10
%   Graphique de derniere trajectoire
figure(1)
clf;
plot(real(zi(:,1)),imag(zi(:,1)),'b.-')
%axis equal
xlim([-5 5])
ylim([-5 5])
hold on
xlabel('x');ylabel('y');
title('Trajectoire');
