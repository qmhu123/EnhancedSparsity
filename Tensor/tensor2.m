% Give the index of a full tensor product basis:
%
%   \phi_a = \phi_a1 \phi_a2 \phi_a3 ... \phi_ad,
%
% where \sum_i^d ai = order.
% 
% Author : Xiu Yang
% Update : July 5th, 2012

function [matrix] = tensor2(dim, order)
matrix = [];
if (dim == 1)
  matrix = order;
else
  if (order == 0)
    matrix = zeros(1, dim);
  else
    for p = 0:order
      sub_mat = tensor2(dim-1, order-p);
      matrix = [matrix; [sub_mat p*ones(size(sub_mat,1), 1)]];
    end
  end
end
