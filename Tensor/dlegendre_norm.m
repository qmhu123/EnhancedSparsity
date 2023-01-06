% Value of the derivative of Legendre polynomial at point x, which can be
% either a scalar or a column vector.
%
% Author: Xiu Yang
% Update: July 5th, 2012

function [dval] = dlegendre_norm(x, order)

val = legendre(x, order+1);
dval = zeros(length(x), order+1);
dval(:,2) = ones(length(x),1);

for i = 2:order
  dval(:,i+1) = (x.*val(:,i+1)-val(:,i))*i./(x.*x-1);
end

ind = find(abs(x-1) < 1e-15);
for j = 2:order
  dval(ind, j+1) = j*(j+1)/2.0*ones(length(ind),1);
end

ind = find(abs(x+1) < 1e-15);
coef = -1;
for j = 2:order
  dval(ind, j+1) = coef*j*(j+1)/2.0*ones(length(ind),1);
  coef = -coef;
end

for i = 1:order
  dval(:, i+1) = dval(:, i+1)*sqrt(2*i+1);
end
