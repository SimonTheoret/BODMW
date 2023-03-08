function matriceQ = gramSchmidtClassiqueItere(matriceInitiale, nbrIterations)
% Adapté pour matrice R^(mxn), n >= m et de rang n

% Itération initiale
if nbrIterations > 1
    conteneur = cell(nbrIterations, 2);
    [Qinit,Rinit] = gramSchmidtClassique(matriceInitiale);
    conteneur{1,1} = Qinit;
    conteneur{1,2} = Rinit;
    
    for i=2:nbrIterations
        Qactuelle = conteneur{i-1,1}; 
        [Qprochaine,Rprochaine]= gramSchmidtClassique(Qactuelle);
        conteneur{i,1} = Qprochaine;
        conteneur{i,2} = Rprochaine;
    end
    matriceQ = Qprochaine;
else
    [matriceQ,pasinteressant] = gramSchmidtClassique(matriceInitiale);
end % modifié si on veut autre chose que la matrice Q