% Value of normalized Charlier polynomial at point x, which can be either a 
% scalar or a column vector.
%
% Author: Xiu Yang
% Update: May 3rd, 2018

function [val] = charlier_norm(x, order, lambda)

val = charlier(x, order, lambda);

if order > 0
  for i = 1:order
    val(:,i+1) = val(:,i+1)/sqrt(lambda^i*factorial(i));
  end
end

