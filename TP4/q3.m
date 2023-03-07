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
% MÉTHODE C - GRAM-SCHMIDT MODIFIÉ
% MÉTHODE D - GRAM-SCHMIDT CLASSIQUE ITÉRÉE
% MÉTHODE E - ROTATIONS DE GIVENS
% MÉTHODE F - TRANSFORMATIONS DE HOUSEHOLDER

% Visualisation sur la console
errA' % A - Cholesky
% B - GM Classique
% C - GM Modifié
% D - GM Classique itérée
% E - Rotation de Givens
% F - Transformations de Householder

%% FIGURE
log10A = -log10(errA);  % A - Cholesky
% B - GM Classique
% C - GM Modifié
% D - GM Classique itérée
% E - Rotation de Givens
% F - Transformations de Householder

figure(1);
hold on;
plot(tailleCho,log10A,'bo--')  % A - Cholesky
% B - GM Classique
% C - GM Modifié
% D - GM Classique itérée
% E - Rotation de Givens
% F - Transformations de Householder

legend('A - Cholesky') % ------------AJOUTER B à F !!!!!!!!!!!!!!!!!!!!!
hold off;
