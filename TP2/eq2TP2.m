function  xpp = eq2TP2(x)
    mu = 0.012155092;
    y=0;
    yp=0;
    xpp =x + 2*yp - (1-mu)*(x+mu)./((x+mu).^2+y.^2).^(3/2) - mu*(x-1+mu)./((x-1+mu).^2+y.^2).^(3/2);

end
