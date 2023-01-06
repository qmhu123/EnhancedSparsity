% Tensor product rule for quadrature points.

% Author : Xiu Yang
% Update : 12/28/2014

function [pnt, wgt] = tensor_rule(dim, pnt1D, wgt1D)

% Number of points in each dimension
num_pnt1D = length(wgt1D);
if (num_pnt1D ~= length(pnt1D))
  display('The number of quadrature points and the weights should be the same');
  return
end

indx_mat = tensor_indx(dim, num_pnt1D);

num_pnt = size(indx_mat,1);
pnt = zeros(num_pnt,dim);
wgt = ones(num_pnt,1);

for k = 1:num_pnt
  for d = 1:dim
    indx = indx_mat(k,d);
    pnt(k,d) = pnt1D(indx);
    wgt(k)   = wgt(k)*wgt1D(indx);
  end
end
