function z = solexFoucault(t,z0,z0Prime,omega,Omega,theta)
sinTheta = sin(theta);
omegaZero = sqrt(omega^2+Omega^2*sinTheta^2);

z = exp(-1i*Omega*sinTheta*t)*(z0(cos(omegaZero*t)+1i*(Omega*sinTheta)/omegaZero*sin(omegaZero*t))+z0Prime/omegaZero*sin(omegaZero*t));
end