function [tt,path]=predcor(odefun,tspan,y0,Nh,predictor,corrector,omega,Omega,theta)
% INPUT
%   odefun    : La fonction de l'ode  
%   tspan     : L'intervalle de temps � r�soudre l'ode
%   y0        : Z�ro pas (Contion initiale)
%   Nh        : Nombre total de pas
%   predictor : Fonction pr�dicteur
%   corrector : Fonction correcteur
% OUTPUT
%   tt        : Intervalle de discr�tisation du temps
%   path      : Trajectoire finale
% -------------------------------------------------------------------------

% Pas de temps h
h=(tspan(2)-tspan(1))/Nh;

% Intervalle de temps selon le nombre de pas
tt=linspace(tspan(1),tspan(2),Nh+1);

% Initialisation des conditions initiales
path=[y0];

% Deux premiers pas avec RK4
[~,uInter]=rk4(odefun,[tspan(1),tspan(1)+2*h],y0,2);
path = [path; uInter(2,:)];
path = [path; uInter(3,:)];
pos=path(end,:);               % Derni�re position

% M�thode pr�dicteur correcteur
for t=tt(3:end-1)
   % Pr�dicteur 
   posPredicteur = feval(predictor, odefun, t, pos, path(end-1,:), path(end - 2, :), h,omega,Omega,theta);
   % Correcteur
   pos = feval(corrector, odefun, t, pos, path(end-1, :), path(end-2, :), posPredicteur, h,omega,Omega,theta);
   
   % Sauvegarde de la nouvelle valeur
   path = [path; pos];
end

end
