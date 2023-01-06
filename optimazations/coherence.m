function [mu,BtB] = coherence(A)
M = size(A,1);
B = A./repmat(sqrt(sum(conj(A).*A)),M,1);
BtB = abs(B'*B);
mu = max(max(BtB-diag(diag(BtB))));
end
