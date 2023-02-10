function [t,u]=predcor(odefun,tspan,y0,y1,y2,Nh,predictor,corrector)
% Il nous faut ajouter y1 et y2 à la matrice u pour appliquer AB3. 

% Taille des pas
h=(tspan(2)-tspan(1))/Nh;

% Initialisation des conditions initiales

w=y2;               % Dernière position
u=[y0,y1,y2].';

% Intervalle de temps selon le nombre de pas
tt=linspace(tspan(1),tspan(2),Nh+1);

% Méthode prédicteur correcteur
for t=tt(3:end-1)
   % Prédicteur 
   upre = feval(predictor, odefun, t, w, u(end-1,:), u(end - 2, :), h);
   % Correcteur
   w = feval(corrector, odefun, t, w, u(end-1, :), u(end-2, :), upre, h);
   
   % Sauvegarde de la nouvelle valeur
   u = [u; w.'];
end

% Intervalle de temps
t = tt;
end
