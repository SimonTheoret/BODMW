function position = AM4Generale (odefun, temps, positionActuelle, positionPrecedente, ...
        avantDernierePosition, positionPredicteur, taillePas)

    function valeurPredicteur = functionPredicteur ()
        valeurPredicteur = feval(odefun, temps + taillePas, positionPredicteur);
    end

    function valeur0 = function0()
        valeur0 = feval(odefun, temps, positionActuelle);
    end

    function valeur1 = function1()
        valeur1 = feval(odefun, temps - taillePas, positionPrecedente);
    end

    function valeur2 = function2()
        valeur2 = feval(odefun, temps - 2 * taillePas, avantDernierePosition);
    end

    function augmentation = increment()
        augmentation = taillePas / 24 * ( 9 * functionPredicteur() + ...
        19 * function0() - 5 * function1() + function2());
    end

    position = positionActuelle + increment();
end


