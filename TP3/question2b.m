%SCRIPT : Question 2 b) du TP3. 
format compact;
format short e;

% Conditions initiales
t0=0;

g = 9.81;
theta=0;
l=67;
Omega = 7.29*10^(-5);
omega = sqrt(g/l);
omegaZero = sqrt(omega^2+Omega^2*sin(theta)^2);

x0=6.0;
x0Prime = 0;
y0 = 0;
y0Prime = 0;
z0 = x0 + 1i * y0;
z0Prime = x0Prime + 1i * y0Prime;
z0Total = [z0,z0Prime];

%Temps finale
tf = 6*pi/omegaZero;
%tf=3*pi/omegaZero;

% Solution de référence
zRef=solexFoucault(tf,z0,z0Prime,omega,Omega,theta)

% Nombre de pas
n=2.^(5:13);
taille = length(n);

% Initialisation de vecteur
zf = zeros(length(taille));
err_pc = zeros(1,length(taille));

% Méthode prédicteur correcteur
for i=1:taille;
ni=n(i); % Nombre de pas nécessaire pour cet itération

% Trajectoire
[ti,zi]=predcor(@foucaultODE,[t0,tf],z0Total,ni,'AB3','AM4',omega,Omega,theta); 
%[ti,zi]=predcor2(@foucaultODE,[t0,tf],z0Total,ni,'eulerunpas2','cnunpas2',omega,Omega,theta); 

zf(i)=zi(end,1);
err_pc(i)=norm((zf(i)-zRef));   % Erreur absolue;
end;
zf
err_pc

%   Graphique
h=(tf-t0)./n;
figure(1)
clf;
polyfit(log(h(1:end)),log(err_pc(1:end)),1)
loglog(h,err_pc,'bo-')
hold on
legend('predicteur-correcteur')
xlabel('h');ylabel('Erreur en norme 2 en t_f');
title("Comportement norme de l'erreur(h) pour predicteur-correcteur AB3 et AM4");

%   Graphique de derniere trajectoire
figure(2)
clf;
plot(real(zi(:,1)),imag(zi(:,1)),'b.-')
hold on
xlabel('x');ylabel('y');
title('Trajectoire');

polyfit(log(h(1:end)),log(err_pc(1:end)),1)
