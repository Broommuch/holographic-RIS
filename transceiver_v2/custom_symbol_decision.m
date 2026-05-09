function dec_symbols = custom_symbol_decision(rx_symbols)
%CUSTOM_SYMBOL_DECISION Classify complex symbols by manually defined regions.
% 调试符号检测用，不用放在主函数里
%
%   dec_symbols = custom_symbol_decision(rx_symbols)
%
%   Input:
%       rx_symbols  - Complex input symbols, can be vector or matrix
%
%   Output:
%       dec_symbols - Detected complex symbols with the same size as rx_symbols

    % Initialize output
    dec_symbols = zeros(size(rx_symbols));

    % Extract real and imaginary parts
    I = real(rx_symbols);
    Q = imag(rx_symbols);

    % Region 1:
    % real > 0.8, imag > 0.2  ->  0.7071 + 0.7071i
    idx0 = (I > 0.8) & (Q > 0.2);
    dec_symbols(idx0) = 0.7071 + 0.7071i;

    % Region 2:
    % 0.4 <= real <= 0.8, 0.2 <= imag <= 0.4  ->  0.7071 - 0.7071i
    idx1 = (I >= 0.4) & (I <= 0.8) & ...
           (Q >= 0.2) & (Q <= 0.4);
    dec_symbols(idx1) = 0.7071 - 0.7071i;

    % Region 3:
    % -0.2 <= real <= 0.4, -0.2 <= imag <= 0.2  ->  -0.7071 + 0.7071i
    idx2 = (I >= -0.2) & (I <= 0.4) & ...
           (Q >= -0.2) & (Q <= 0.2);
    dec_symbols(idx2) = -0.7071 + 0.7071i;

    % Other regions:
    % remaining points -> -0.7071 - 0.7071i
    idx_other = ~(idx0 | idx1 | idx2);
    dec_symbols(idx_other) = -0.7071 - 0.7071i;

end