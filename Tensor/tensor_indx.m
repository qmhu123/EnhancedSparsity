% Generate the matrix of the indices for tensor product rule

% Author : Xiu Yang
% Update : 12/28/2014

function mat = tensor_indx(dim, num_pnt)
mat = zeros(num_pnt^dim, dim);
if (dim == 1)
  mat = (1:num_pnt)';
else
  sub_mat = tensor_indx(dim-1, num_pnt);
  row = size(sub_mat,1);
  for k = 1:num_pnt
    mat((k-1)*row+1:k*row,:) = [k*ones(row,1) sub_mat];
  end
end
