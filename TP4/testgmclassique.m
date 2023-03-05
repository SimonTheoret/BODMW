% test GMclassique. Ça morche
A = [3,2; 1,2];
[Q,R]=gramSchmidtClassique(A);
%Q*R
disp("test itéré:")
% test GM itéré
sol = gramSchmidtClassiqueItere(A,1);
autSol = gramSchmidtClassiqueItere(A,7)

% gwo test de GM itéré:
