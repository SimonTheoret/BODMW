% Exemple en classe de prédicteur (Euler) correcteur (CN)
clear;clc;
format compact;
format short e;
% Intervalle de temps
t0=0;tf=20;

% Condition initiale
y0=0.5;

% Solution de référence
yref=solex(tf);

% Nombre de pas
taille = 1:15;
n=2.^taille;

% Initialisation de vecteur
yf = zeros(1,length(taille));
err_pc = zeros(1,length(taille));

% Méthode prédicteur correcteur
for i=taille;
ni=n(i); % Nombre de pas nécessaire pour cet itération


% Deux premiers pas avec RK4
h=(tf-t0)/ni;
[tInter,uInter]=rk4(@f1,[t0,t0+2*h],y0,2);
y1 = uInter(2,:);
y2 = uInter(3,:);
% Trajectoire
[ti,yi]=predcor(@f1,[t0,tf],y0,y1,y2,ni,'AB3','AM4'); 

yf(i)=yi(ni+1);
err_pc(i)=abs((yf(i)-yref)/yref);   % Erreur relative;
end;

%   Graphique
h=(tf-t0)./n;
figure(1)
polyfit(log(h(5:end)),log(err_pc(5:end)),1)
loglog(h,err_pc,'bo-')
hold on
legend('predicteur-correcteur')
xlabel('h');ylabel('Erreur relative en t_f');
title('Comportement erreur(h) pour predicteur-correcteur Euler explicite-Crank Nicolson');