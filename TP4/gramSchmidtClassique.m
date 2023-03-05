function [matriceQ, matriceR] = gramSchmidtClassique(matriceInitiale)
%Adapté pour matrice R^(mxn), n >= m et de rang n

rang = size(matriceInitiale,2);
tailleMatriceQ = size(matriceInitiale);
taileMatriceR = rang;

matriceQ = zeros(tailleMatriceQ);
matriceR = zeros(taileMatriceR);

for j = 1:rang
    vj = matriceInitiale(:,j);
    for i = 1:j-1
        matriceR(i,j) = dot(matriceQ(:,i), matriceInitiale(:,j));
        vj = vj - matriceR(i,j) * matriceQ(:,i);
    end
    matriceR(j,j) = norm(vj, 2);
    matriceQ(:,j) = vj / matriceR(j,j);
end
