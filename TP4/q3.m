% Question 3
clc; clear;
format compact;
format short e;
%Dimensions du problèmes
m=50;   
taille = 4:18;
p=2;

% MÉTHODE A - CHOLESKY
% Initialisation vecteur erreur
tailleCho = 4:12;
errA = zeros(length(tailleCho),1);
for n=tailleCho;
    % Générer le vecteur des ti et la matrice A associé
    [t,A] = genererTiMatrice(m,n);

    % Définir le xExact qu'on veut retrouver
    xExact = ones(1,n);

    % Former notre b par la multiplication
    b = A*xExact';

    % Cholesky pour retrouver x
    ACho = (A')*A;
    bCho = A'*b;
    L = cholesky(ACho);
    qChap = A*inv(L');
    
    errA(n-tailleCho(1)+1)=norm(eye(n) - qChap'*qChap);
%     z = subAvant(L,bCho);
%     xCho = subArriere(L',z);
end

% MÉTHODE B - GRAM-SCHMIDT CLASSIQUE
errB = zeros(length(taille), 1);
for n=taille
    [t,A] = genererTiMatrice(m,n);

    [qChap,R]=gramSchmidtClassique(A);

    errB(n-taille(1)+1)=norm(eye(n) - qChap'*qChap);
end
% MÉTHODE C - GRAM-SCHMIDT MODIFIÉ
errC = zeros(length(taille),1);
for n=taille
    % Générer le vecteur des ti et la matrice A associé
    [t,A] = genererTiMatrice(m,n);

    [qChap,R] = mgs(A);
    
    errC(n-taille(1)+1)=norm(eye(n) - qChap'*qChap);

end
% MÉTHODE D - GRAM-SCHMIDT CLASSIQUE ITÉRÉE
errD = zeros(length(taille), 1);
for n=taille
    [t,A] = genererTiMatrice(m,n);

    qChap=gramSchmidtClassiqueItere(A,2);

    errD(n-taille(1)+1)=norm(eye(n) - qChap'*qChap);
end
% MÉTHODE E - ROTATIONS DE GIVENS
errE = zeros(length(taille),1);
for n=taille
    % Générer le vecteur des ti et la matrice A associé
    [t,A] = genererTiMatrice(m,n);

    [qChap,R] = rotationsGivens(A);
    
    errE(n-taille(1)+1)=norm(eye(n) - qChap'*qChap);

end
% MÉTHODE F - TRANSFORMATIONS DE HOUSEHOLDER
errF = zeros(length(taille),1);
for n=taille
    % Générer le vecteur des ti et la matrice A associé
    [t,A] = genererTiMatrice(m,n);

    [qChap,R] = transfHouseholder(A);
    eye(n) - qChap'*qChap
    errF(n-taille(1)+1)=norm(eye(n) - qChap'*qChap);

end

% Visualisation sur la console
errA' % A - Cholesky

% B - GM Classique
errB'

% C - GM Modifié
errC'

% D - GM Classique itérée
errD'

% E - Rotation de Givens
errE' 

% F - Transformations de Householder
errF'

%% FIGURE
% A - Cholesky
log10A = -log10(errA);  

% B - GM Classique
log10B = -log10(errB);

% C - GM Modifié
log10C = -log10(errC);

% D - GM Classique itérée
log10D = -log10(errD);

% E - Rotation de Givens
log10E = -log10(errE); 

% F - Transformations de Householder
log10F = -log10(errF); 

figure(1);
hold on;
% A - Cholesky
plot(tailleCho,log10A,'bo--')  

% B - GM Classique
plot(taille,log10B,'gx--')  

% C - GM Modifié
plot(taille,log10C,'ko--') 

% D - GM Classique itérée
plot(taille,log10D,'cx-') 

% E - Rotation de Givens
plot(taille,log10E,'r*--') 

% F - Transformations de Householder
%plot(taille,log10E,'r*--') 

legend('A - Cholesky','C - GM Modifié','E - Rotation de Givens','Location','best') % ------------AJOUTER B à F !!!!!!!!!!!!!!!!!!!!!
hold off;
