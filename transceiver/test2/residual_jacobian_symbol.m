function [r,J] = residual_jacobian_symbol(x,V,A,b,I)
L = length(x)/2;
s = x(1:L) + 1j*x(L+1:end);
y = V*A*s + b;
mu = abs(y).^2;
r = mu - I;

eps = 1e-6;
K = length(I);
J = zeros(K,2*L);

for k = 1:2*L
    xp = x; xp(k) = xp(k) + eps;
    sp = xp(1:L) + 1j*xp(L+1:end);
    yp = V*A*sp + b;
    mup = abs(yp).^2;
    J(:,k) = (mup - mu)/eps;
end
end