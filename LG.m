function U = LG(R,Phi,P,L,weights,w0)
% This function computes a superposition of LG mode at the plane z=0
% R,Phi are coordinate matrices
% "weights" is a weight vector for the coefficients in the superposition
% E.g. generate mode U = 5*LG_0^-5 + I*LG_0^3 - LG_1^4 + 2*LG_3^4 then,
% x = -1:0.01:1, y = -1:0.01:1, [X,Y] = meshgrid(x,y); [Phi,R] = cart2pol(X,Y);
% P = [0,0,1,3], L = [-5,3,4,4], weights = [5,I,-1,2]; w0 = 1;
% U = LG(R,Phi,P,L,weights,w0);

U = zeros(size(R)); % initialise field

for i = 1:length(weights)
    U = U + weights(i).*(sqrt(2*factorial(P(i))/(pi*factorial(P(i)+abs(L(i))))).*(1/w0).*(sqrt(2).*R./w0).^(abs(L(i))) .* exp(-R.^2./w0^2).*Laguerre(P(i),abs(L(i)),2.*R.^2./(w0.^2)).* exp(1i.*L(i).*Phi));
end

end

function y = Laguerre(p,l,x)
% computes associated Laguerre polynomials 
if p == 0
    y = ones(size(x));
elseif p == 1
    y = 1 + l - x;
else
    y = ((2*p+l-1-x)./p).*Laguerre(p-1,l,x) - ((p+l-1)./p).*Laguerre(p-2,l,x);
end
end
