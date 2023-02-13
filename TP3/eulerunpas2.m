function [u] = eulerunpas2(t,y,h,f)
u = y + h*f;
return
