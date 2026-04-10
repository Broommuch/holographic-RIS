function s_hat = estimate_symbol_GN(V,A,b,I)
L = size(V,2);
x = zeros(2*L,1);
for iter = 1:20
    [r,J] = residual_jacobian_symbol(x,V,A,b,I);
%     dx = -(J.'*J)\(J.'*r);
    mu = 1e-3 * trace(J.'*J) / size(J,2);   % 自适应阻尼
    dx = -(J.'*J + mu*eye(size(J,2))) \ (J.'*r);
    x = x + dx;
    if norm(dx) < 1e-6, break; end
end
s_hat = x(1:L) + 1j*x(L+1:end);
end