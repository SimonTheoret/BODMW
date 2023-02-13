function [tt,u]=rk4TP3(odefun,tspan,y0,Nh)
%
% script MATLAB écrit par Anne Bourlioux
%
tt=linspace(tspan(1),tspan(2),Nh+1);

h=(tspan(2)-tspan(1))/Nh;  hh=h/6;

u=y0; %ex: y0 = [x1,x2,...,xn]
      %y0 est un vecteur ligne de valeurs init2iales

for t=tt(1:end-1)
  y = u(end,:); %toute la dernière ligne
  k1=feval(odefun,t,y);
  y1 = y + 0.5*h*k1;
  t1 = t+h/2;
  k2=feval(odefun,t1,y1);
  y2 = y + 0.5*h*k2;
  t2 = t+h/2;
  k3=feval(odefun,t2,y2);
  y3 = y + h*k3;
  t3=t+h;
  k4=feval(odefun,t3,y3);
  %
  u = [u; u(end,:) + hh*(k1+2*k2+2*k3+k4)]; %Matrice Nhx4
end