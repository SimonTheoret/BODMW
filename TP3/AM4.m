function position = AM4 (odefun, temps, positionActuelle, positionPrecedente, ...
        avantDernirePosition, positionPredicteur, taillePas)

    funcPredicteur = feval(odefun, temps + taillePas, positionPredicteur);

    func0 = feval(odefun, temps, positionActuelle);

    func1 = feval(odefun, temps - taillepas, positionPrecedente);

    func2 = feval(odefun, temps - 2 * taillePas, avantDernirePosition);

    increment = taillePas / 24 * ( 9 * funcPredicteur + 19 * func0 - ...
        5 * func1 + func2);

    position = positionActuelle + increment;
end
