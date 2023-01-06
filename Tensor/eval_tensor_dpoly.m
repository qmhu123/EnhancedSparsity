% Evaluate the derivative polynomial of tensor product form.
%
% Input :
%     poly  : A function of polynomial returns a vector with the values of the
%             polynomials at x according to the index matrix ind_mat.
%     x     : Point to evaluate the polynomial.
%     order : Order of the polynomial.
%     ind_mat: Index matrix.
%
% Output :
%     dval  : value of the derivative of the polynomial, which is a matrix
%
% Author : Xiu Yang
% Update : July 5th, 2012

function [dval] = eval_tensor_dpoly(polynomial, dpolynomial, x, dim, order, ind_mat)

[row, col] = size(ind_mat);
dval = zeros(row, col);

val_mat = polynomial(x, order);
dval_mat = dpolynomial(x, order);

for i = 1:row
  indices = find(ind_mat(i,:));
  for j = 1:length(indices)
    tmp = dval_mat(indices(j), ind_mat(i, indices(j))+1);
    for k = 1:j-1
      tmp = tmp*val_mat(indices(k), ind_mat(i, indices(k))+1);
    end
    for k = j+1:length(indices)
      tmp = tmp*val_mat(indices(k), ind_mat(i, indices(k))+1);
    end
    dval(i, indices(j)) = tmp;
  end
end
