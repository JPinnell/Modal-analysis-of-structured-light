function U = HG(X,Y,N,M,weights,w0)
% This function computes a superposition of HG mode at the plane z=0
% X,Y are matrices, e.g. x = -1:0.1:1, y = -1:0.1:1, [X,Y] = meshgrid(x,y)
% "weights" is a weight vector for the coefficients in the superposition
% E.g. generate mode U = 5*HG_0^-5 + I*HG_0^3 - HG_1^4 + 2*HG_3^4 then,
% p = [0,0,1,3], l = [-5,3,4,4], weights = [5,I,-1,2]

U = zeros(size(X)); %initialise electric field
 
for i = 1:length(weights)
    U = U + weights(i).*(1./w0).*(sqrt(2.^(1-N(i)-M(i))./(pi.*factorial(N(i)).*factorial(M(i)))).*Hermite(N(i),sqrt(2).*X./w0).*Hermite(M(i),sqrt(2).*Y./w0).*exp(-(X.^2+Y.^2)./w0.^2));
end

end

function y = Hermite(n,x)
% Computes Hermite polynomials
if n == 0
    y = ones(size(x));
elseif n == 1
    y = 2.*x;
else
    y = 2.*x.*Hermite(n-1,x) - 2.*(n-1).*Hermite(n-2,x);
end

end