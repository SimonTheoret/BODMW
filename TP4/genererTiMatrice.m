function [t,A] = genererTiMatrice(m,n)
t = zeros(1,m);
A = ones(m,n);
for i=1:m
    t(i) = (i-1)/(m-1);
    for j=2:n
        A(i,j) = t(i)^(j-1);
    end
end
end