function zf = fTroisCorps(z)
% DESCRITION
%       Partie droite du système d'équations
% INPUT
%       z       : Vecteur colonne(4x1) contenant les anciennes valeurs
%
% OUTPUT
%       zf      : Vecteur final [x,y,x',y']
%
% ---------------------------------------------------------------------

mu = 0.012155092;
zf = zeros(4,1);
x = z(1);
y = z(2);
xp = z(3);
yp = z(4);

zf(1) = xp;
zf(2) = yp;
zf(3) = x + 2*yp - (1-mu)*(x+mu)/((x+mu)^2+y^2)^(3/2) - mu*(x-1+mu)/((x-1+mu)^2+y^2)^(3/2);
zf(4) = y - 2*xp - (1-mu)*y/((x+mu)^2+y^2)^(3/2) - mu*y/((x-1+mu)^2+y^2)^(3/2);

end

