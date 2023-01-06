% Value of the derivative of Charlier polynomial at point x, which can be
% either a scalar or a column vector.
%
% Author: Xiu Yang
% Update: May 3rd, 2018

function [dval] = dcharlier_norm(x, order, lambda)

dval = dcharlier(x, order, lambda);
 
if order > 0
  for i = 1:order
    dval(:,i+1) = dval(:,i+1)/sqrt(lambda^i*factorial(i));
  end
end
