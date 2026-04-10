function [r,J] = residual_jacobian_doa(x,I,S,b,Ny,Nz,dy,dz,lambda)

L = size(S,1); T = size(S,2);
[theta,phi,alpha] = unpack_param(x,L);
K = Ny*Nz;

r = zeros(K*T,1);
J = zeros(K*T,length(x));
idx = 1;

for t = 1:T
    y = zeros(K,1);
    for l = 1:L
        v = steering2D(theta(l),phi(l),Ny,Nz,dy,dz,lambda);
        y = y + alpha(l)*v*S(l,t);
    end
    e = abs(y+b).^2 - I(:,t);
    r(idx:idx+K-1) = e;

    eps = 1e-6;
    for p = 1:length(x)
        xp = x; xp(p) = xp(p) + eps;
        yp = zeros(K,1);
        [thp,php,alp] = unpack_param(xp,L);
        for l = 1:L
            vp = steering2D(thp(l),php(l),Ny,Nz,dy,dz,lambda);
            yp = yp + alp(l)*vp*S(l,t);
        end
        ep = abs(yp+b).^2 - I(:,t);
        J(idx:idx+K-1,p) = (ep - e)/eps;
    end
    idx = idx + K;
end
end