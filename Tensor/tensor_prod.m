% Provide the index of tensor product rule
%
% Author : Xiu Yang
% Update : July 5th, 2012

function [matrix] = tensor_prod(dim, n)
matrix = [];

if (n < 1)
  display('Number of index in each dimension should be larger than 0!');
elseif (n == 1)
  matrix = ones(1,dim);
else 
  if (dim < 1)
    display('Dimension should be larger than 0!');
  elseif (dim == 1)
    matrix = (1:n)';
  else
    tmp = tensor_prod(dim-1, n);
    [row, col] = size(tmp);
    for k = 1:n
      matrix = [matrix; [tmp ones(row,1)*k]];
    end
  end
end
