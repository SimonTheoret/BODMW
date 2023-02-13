function [t,u]=predcor2(odefun,tspan,y0,Nh,predictor,corrector,omega,Omega,theta)
%
% OUTPUT
% u : Vecteur des trajectoires
% t : Vecteur des temps

% Taille de pas
h=(tspan(2)-tspan(1))/Nh;

% Initialisation des conditions initiales
y=y0(:); 
w=y;
u=y.';

tt=linspace(tspan(1),tspan(2),Nh+1);

for t=tt(1:end-1)
   fk = feval(odefun,t,w,omega,Omega,theta);
   upre = feval(predictor,t,w,h,fk);
   w = feval(corrector,t+h,w,upre,h,odefun,fk,omega,Omega,theta);
   u = [u; w.'];
end

t = tt;
end
