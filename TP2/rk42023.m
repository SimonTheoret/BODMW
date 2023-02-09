function [Z,tList] = rk42023(ff,intervalle,z0,N)
% INPUT
%    ff      : Fonction qui décrit le système d'EDO transformé;
% intervalle : Intervalle de temps de la simulation;
%    z0      : Position initiale du satellite;
%    N       : Nombre de pas de temps. 
%
% OUTPUT
%    Z       :  Tajectoire du satellite (Z = [x,y,x',y']) (Vecteur (4x1));
%    tList   :  Liste des temps t0,t1,...tf.
%
% ---------------------------------------------------------------------
% INITIALISATION

% Temps initiale et final
t0 = intervalle(1);
tf = intervalle(2);
    
% Pas de temps initial selon Fortin743
h = (tf-t0)/N;

% Initialiation condition initiale
Z = zeros(4,N+1);
Z(:,1) = z0;

tList = zeros(1,N+1);
tList(1) = t0;

% ---------------------------------------------------------------------
% MÉTHODE NUMÉRIQUE RUNGE-KUTTA ORDRE 4 (RK4)
for i=1:N
    
    path = Z(:,i);         % État du système précédent au temps précédent
    
    % MÉTHODE RUNGE-KUTTA D'ORDRE 4
    p1 = path;
    k1 = h*feval(ff,p1);
    p2 =path+1/4*k1;
    k2= h*feval(ff,p2);
    p3 = path+3/32*k1+9/32*k2;
    k3= h*feval(ff,p3);
    p4 = path+1932/2197*k1-7200/2197*k2+7296/2197*k3;
    k4= h*feval(ff,p4);
    p5 = path+439/216*k1-8*k2+3680/513*k3-845/4104*k4;
    k5= h*feval(ff,p5);

 
    Z(:,i+1)  = path+(25/216*k1+1408/2565*k3 +2197/4104*k4-1/5*k5);
    
    % SAUVEGARDE DE DONNÉES AVANT LE PROCHAIN PAS DE TEMPS
    tList(i+1) = tList(i)+h;   % Sauvegarde de ce temps
end

end

