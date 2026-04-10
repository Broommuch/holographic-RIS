function s = qpsk_hard(s)
s = sign(real(s)) + 1j*sign(imag(s));
s(real(s)==0) = 1 + 1j*sign(imag(s(real(s)==0)));
s(imag(s)==0) = sign(real(s(imag(s)==0))) + 1j;
end