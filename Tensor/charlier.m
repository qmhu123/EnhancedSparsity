% Value of Charlier polynomial at point x, which can be either a scalar or a
% column vector.
%
% Author: Xiu Yang
% Update: May 3rd, 2018

function [val] = charlier(x, order, lambda)

val = ones(length(x), order+1);

if order > 0
  for i = 1:order
    val(:,i+1) = charlierF(x, i, lambda);
  end
end

