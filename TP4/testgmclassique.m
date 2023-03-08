% test GMclassique. Ça morche
A = [3,2; 1,2];
%[Q,R]=gramSchmidtClassique(A);
%Q*R
%disp("test itéré:")
% test GM itéré
%sol = gramSchmidtClassiqueItere(A,1);
%autSol = gramSchmidtClassiqueItere(A,7);

% test transformation householder:
%transfHouseholder(A);

B = [12,-51,4; 6, 167, -68; -4, 24 ,-41];
disp("Boy HH:")
[Q,R]=transfHouseholder(B)
Q*R
%transfHouseholder(A)

%C = [4,1,-2,2; 1,2,0,1; -2,0,3,-2; 2,1,-2,-1];
%[Q,R]=transfHouseholder(C)
