function [tt,path]=predCorrGenerale(odefun,tspan,y0,Nh,predictor,corrector)
% INPUT
%   odefun    : La fonction de l'ode  
%   tspan     : L'intervalle de temps à résoudre l'ode
%   y0        : Zéro pas (Contion initiale)
%   Nh        : Nombre total de pas
%   predictor : Fonction prédicteur
%   corrector : Fonction correcteur
% OUTPUT
%   tt        : Intervalle de discrétisation du temps
%   path      : Trajectoire finale
% -------------------------------------------------------------------------

% Pas de temps h
h=(tspan(2)-tspan(1))/Nh;

% Intervalle de temps selon le nombre de pas
tt=linspace(tspan(1),tspan(2),Nh+1);

% Initialisation des conditions initiales
path=[y0];

% Deux premiers pas avec RK4
[~,uInter]=rk4TP3(odefun,[tspan(1),tspan(1)+2*h],y0,2,omega,Omega,theta);
path = [path; uInter(2,:)];
path = [path; uInter(3,:)];
pos=path(end,:);               % Dernière position

% Méthode prédicteur correcteur
for t=tt(3:end-1)
   % Prédicteur 
   posPredicteur = feval(predictor, odefun, t, pos, path(end-1,:), path(end - 2, :), h);
   % Correcteur
   pos = feval(corrector, odefun, t, pos, path(end-1, :), path(end-2, :), posPredicteur, h);
   
   % Sauvegarde de la nouvelle valeur
   path = [path; pos];
end

end
