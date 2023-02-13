function z = solexFoucault(t,z0,z0Prime,omega,Omega,theta)
sinTheta = sin(theta);
omegaZero = sqrt(omega^2+Omega^2*sinTheta^2);
p1 = exp(-1i*Omega*sinTheta*t);
p2 = cos(omegaZero*t)+1i*(Omega*sinTheta)/omegaZero*sin(omegaZero*t);
p3 = z0Prime/omegaZero*sin(omegaZero*t);
z = p1*(z0*p2+p3);
end