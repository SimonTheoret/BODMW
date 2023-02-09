function [Z,hList,tList] = rkf2023(ff,intervalle,z0,tol)
% INPUT
%    ff      : Fonction qui décrit le système d'EDO transformé;
% intervalle : Intervalle de temps de la simulation;
%    z0      : Position initiale du satellite;
%    tol     : Tolérance choisie sur la précision. 
%
% OUTPUT
%    Z       :  Tajectoire du satellite (Z = [x,y,x',y']) (Vecteur (4x1));
%    hList   :  Liste des pas de temps h1,h2,...,hf;
%    tList   :  Liste des temps t0,t1,...tf.
%
% ---------------------------------------------------------------------
% INITIALISATION

% Temps initiale et final
t0 = intervalle(1);
tf = intervalle(2);

% Pas de temps initial selon Fortin743
hmax = (tf-t0)/16;
h = hmax/8;

% Initialiation condition initiale
Z(:,1) = z0;
t = t0;
tList(1) = t;

% ---------------------------------------------------------------------
% MÉTHODE NUMÉRIQUE Runge-Kutta-Fehlberg (RKF)
while t < tf
    
    path = Z(:,end);         % État du système précédent au temps précédent
    Z(:,end+1) = zeros(4,1); % Allocation mémoire pour prochain pas
    
    % CALCUL DE L'ERREUR
    % Calcul des constantes
    k1 = h*feval(ff,path);                                  
    p2 = path+(1/4)*k1;
    k2 = h*feval(ff,p2);
    p3 = path+(3/32)*k1+(9/32)*k2;
    k3 = h*feval(ff,p3);
    p4 = path+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3;
    k4 = h*feval(ff,p4);
    p5 = path+(439/216)*k1-(8)*k2+(3680/513)*k3-(845/4104)*k4;
    k5 = h*feval(ff,p5);
    p6 = path-(8/27)*k1+(2)*k2-(3544/2565)*k3+(1859/4104)*k4-(11/40)*k5;
    k6 = h*feval(ff,p6);
    
    % Erreur entre RKF4 et RKF5
    Emat = (1/360)*k1-(128/4275)*k3-(2197/75240)*k4+(1/50)*k5+(2/55)*k6;
    E = norm(Emat.','inf');                  
    
    % Facteur d'augemntation beta
    beta = (tol/2/E)^(1/5); 
    
    % MODIFICATION DU PAS DE TEMPS
    % Vérification de ne pas dépasser le pas de temps maximal
    if (beta*h)>hmax 
        h = hmax;
    else
        h  = beta*h;
    end

    % Vérification pour ne pas dépasser le temps final
    if t+h > tf
        h=tf-t;
    end
    
    
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

 
    Z(:,end)  = path+(25/216*k1+1408/2565*k3 +2197/4104*k4-1/5*k5);
    
    % SAUVEGARDE DE DONNÉES AVANT LE PROCHAIN PAS DE TEMPS
    % Sauvegarde du pas de temps
    if t==t0
        hList(1) = h;
    else
        hList(end+1) = h;       
    end
    % Sauvegarde du temps   
    t = t+h;            
    tList(end+1) = t;   
end

end