function [matriceQ,matriceR] = transfHouseholder(matriceInitiale)
nbrColonnes = size(matriceInitiale,2);
nbrLignes = size(matriceInitiale,1);

matriceQ = eye(nbrLignes,nbrColonnes);

    for k=1:nbrColonnes
        x = matriceInitiale(k:end, k);
        matIdentite = eye(size(x,1));
        vectCanonique = matIdentite(:,1);
        vk = sign(x(1))*norm(x,2)*vectCanonique + x;
        vk = vk/norm(vk,2);
        matriceInitiale(k:end,k:end)= matriceInitiale(k:end,k:end)- 2*vk*(vk'*matriceInitiale(k:end,k:end));
        matriceQPrime = matIdentite - 2*vk*vk';
        matriceQPrime = blkdiag(eye(nbrLignes-size(x,1)),matriceQPrime);
        matriceQ = matriceQ*matriceQPrime;
    end
    matriceR = matriceInitiale;
end