function [tt,u]=rk4(odefun,tspan,y0,Nh)
%
% script MATLAB �crit par Anne Bourlioux
%
tt=linspace(tspan(1),tspan(2),Nh+1);

h=(tspan(2)-tspan(1))/Nh;  hh=h/6;

u=y0; %ex: y0 = [x1,x2,...,xn]
      %y0 est un vecteur ligne de valeurs initiales

for t=tt(1:end-1)
  y = u(end,:); %toute la derni�re ligne
  k1=feval(odefun,y);
  y1 = y + 0.5*h*k1;
  k2=feval(odefun,y1);
  y2 = y + 0.5*h*k2;
  k3=feval(odefun,y2);
  y3 = y + h*k3;
  k4=feval(odefun,y3);
  %
  u = [u; u(end,:) + hh*(k1+2*k2+2*k3+k4)]; %Matrice Nh x 4
end