function positionFinale = foucaultODE(temps,positionInitiale,omega,Omega,theta)
% DESCRITION
%       Partie droite du système d'équations
% INPUT
%       positionInitiale    : Vecteur colonne(2x1) contenant les anciennes
%       valeurs z1,z2
%
% OUTPUT
%       positionFinale      : Vecteur final [z,z'] = [z1,z2]
%
% ---------------------------------------------------------------------
z1= positionInitiale(1);
z2= positionInitiale(2);

positionFinale(1) = z2;

positionFinale(2) = (-omega^2)*z1-(2i*Omega*sin(theta))*z2;

end

