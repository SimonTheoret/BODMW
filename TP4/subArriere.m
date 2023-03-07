function x = subArriere(A,b)

%Dimensions A et dimensions b
[nbLigneA,nbColonneA] = size(A);
nbLigneb = length(b);

% Test si dimensions matrice A et vecteur b concordent
if  nbLigneA ~= nbLigneb
    disp("Dimensions de matrice et vecteur incompatible.")
    return
end

% Initialisation x
x = zeros(nbLigneb,1);

% Substitution arrière
for j=nbLigneA:-1:1
    x(j) = (b(j)-dot(A(j,j+1:end),x(j+1:end)))/A(j,j);
end