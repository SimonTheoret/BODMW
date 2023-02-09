function [t,u]=predcor(odefun,tspan,y0,Nh,predictor,corrector)
%  PREDCOR  R�sout une �quation diff�rentielle avec une
%  m�thode predicteur-correcteur
%  [T,Y]=PREDCOR(ODEFUN,TSPAN,Y0,NH,PRED,CORR) avec
%  TSPAN=[T0 TF]
%  int�gre le syst�me d'�quations diff�rentielles
%  y'=f(t,y) du temps T0 au temps TF avec la condition
%  initiale Y0 en utilisant une m�thode g�n�rale
%  pr�dicteur-correcteur sur une grille de NH
%  intervalles �quidistribu�s. La fonction ODEFUN(T,Y)
%  doit retourner un vecteur correspondant � f(t,y)
%  de m�me dimension que Y.
%  Chaque ligne de la solution Y correspond
%  � un temps du vecteur colonne T.
%  [T,Y]=PREDCOR(ODEFUN,TSPAN,Y0,NH,PRED,CORR,P1,..)
%  passe les param�tres suppl�mentaires P1,P2,.. aux
%  fonctions ODEFUN, PRED et COOR de la mani�re
%  suivante: ODEFUN(T,Y,P1,...), PRED(T,Y,P1,P2...),
%  CORR(T,Y,P1,P2...).

% Il nous faut ajouter y1 et y2 � la matrice u pour appliquer AB3. 

h=(tspan(2)-tspan(1))/Nh;
y=y0(:); w=y; u=y.';
tt=linspace(tspan(1),tspan(2),Nh+1);
for t=tt(1:end-1)
   upre = feval(predictor, odefun, t, w, u(end,:), u(end - 1, :), h);
   w = feval(corrector, odefun, t, w, u(end, :), u(end-1, :), upre, h);
   u = [u; w.'];
end
t = tt;
end
