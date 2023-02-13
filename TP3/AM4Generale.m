function position = AM4 (odefun, temps, positionActuelle, positionPrecedente, ...
        avantDernierePosition, positionPredicteur, taillePas,omega,Omega,theta)

    function valeurPredicteur = functionPredicteur ()
        valeurPredicteur = feval(odefun, temps + taillePas, positionPredicteur,omega,Omega,theta);
    end

    function valeur0 = function0()
        valeur0 = feval(odefun, temps, positionActuelle,omega,Omega,theta);
    end

    function valeur1 = function1()
        valeur1 = feval(odefun, temps - taillePas, positionPrecedente,omega,Omega,theta);
    end

    function valeur2 = function2()
        valeur2 = feval(odefun, temps - 2 * taillePas, avantDernierePosition,omega,Omega,theta);
    end

    function augmentation = increment()
        augmentation = taillePas / 24 * ( 9 * functionPredicteur() + ...
        19 * function0() - 5 * function1() + function2());
    end

    position = positionActuelle + increment();
end
