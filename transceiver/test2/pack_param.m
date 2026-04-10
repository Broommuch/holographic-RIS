function x = pack_param(theta,phi,alpha)
x = [theta(:); phi(:); real(alpha); imag(alpha)];
end