function [matriceR, vectInitialTransforme] = householderSolver(matriceInitiale, vectInitial)
nbrLignes = size(matriceInitiale,1);
nbrColonnes = size(matriceInitiale,2);


    for k=1:nbrColonnes
        x = matriceInitiale(k:nbrLignes, k);
        matIdentite = eye(size(x,1));
        vectCanonique = matIdentite(:,1);
        vk = sign(x(1))*norm(x,2)*vectCanonique + x;
        vk = vk/norm(vk,2);
        matriceInitiale(k:nbrLignes,k:nbrColonnes)= matriceInitiale(k:nbrLignes,k:nbrColonnes)- 2*vk*(vk'*matriceInitiale(k:nbrLignes,k:nbrColonnes));
        vectInitial(k:nbrColonnes) = vectInitial(k:nbrColonnes)-2*vk(vk'*vectInitial(k:nbrColonnes));
    end
    matriceR = matriceInitiale;
    vectInitialTransforme = vectInitial;
end