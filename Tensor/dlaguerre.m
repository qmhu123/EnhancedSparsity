% Value of derivatives of Laguerre polynomial at point x, which can be either
% a scalar or a column vector.
%
% Author: Xiu Yang
% Update: August 5th, 2012

function [dval] = dlaguerre(x, order)

dim = length(x);
dval = zeros(dim, order+1);
val = laguerre(x, order);

non0_ind = find(x~=0);
zero_ind = setdiff(1:dim, non0_ind);
tmp = ones(length(zero_ind),1);

for i = 2:order+1
  dval(non0_ind,i) = (i-1)*(val(non0_ind,i) - val(non0_ind,i-1))./x(non0_ind);
  dval(zero_ind,i) = -(i-1)*tmp;
end

