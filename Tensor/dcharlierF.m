% Value of derivative of Charlier polynomial at point x, which can be either a
% scalar or a column vector.
%
% Author: Xiu Yang
% Update: May 3rd, 2018

function [dval] = dcharlierF(x, order, lambda)

if order == 0
  dval = zeros(length(x), 1);
end

if order == 1
  dval = ones(length(x), 1);
end

if order > 1
  dval = x.*dcharlierF(x-1, order-1, lambda) + charlierF(x-1, order-1, lambda) - lambda*dcharlierF(x, order-1, lambda);
end

