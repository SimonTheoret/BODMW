% Méthode Adams-Bashfort (AB3) comme prédicteur dans un schéma numérique
% prédicteur-correcteur la sortie doit être un vecteur

function [position] = AB3Generale( func, temps, positionActuelle,...
        positionPrecedente, avantDernierePosition, taillePas)

    func0 = feval( func, temps, positionActuelle);

    func1 = feval( func, temps - taillePas, positionPrecedente);

    func2 = feval( func, temps - 2 * taillePas, avantDernierePosition);

    increment = taillePas / 12 * ( 23 * func0 - 16 * func1 + 5 * func2 );

    position = positionActuelle + increment;
end

