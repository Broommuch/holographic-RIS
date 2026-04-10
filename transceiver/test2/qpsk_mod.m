function S = qpsk_mod(L,T)
bits = randi([0 1],2*L,T);
S = (2*bits(1:2:end,:) - 1) ...
  + 1j*(2*bits(2:2:end,:) - 1);
end