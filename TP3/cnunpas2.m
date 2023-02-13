function [u] = cnunpas2(t,u,y,h,f,fn,omega,Omega,theta)
u = u + 0.5*h*(feval(f,t,y,omega,Omega,theta)+fn);
return