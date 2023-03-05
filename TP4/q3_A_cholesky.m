% Question 3 pour l'algo A - Cholesky
m=6;
n=5;
[t,A] = genererTiMatrice(m,n);
xExact = ones(1,n);
b = A*xExact';
ACho = (A')*A;
cholesky(ACho)