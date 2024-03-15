function [M] = vec2matx(V, N)

format long;

V = V.';
V = V(:);

R = ceil(length (V) / N);

M = reshape (V, N, R).';

end %endfunction