clear;
clc;

%f = @(x) exp(x(:,1)).*sin(x(:,2)) + x(:,3).^2;
%f = @(x) exp(x(:,1));
f = @(x) 10.^(x(:,1));

dim = 3;
poly_order = 6;
indx_mat = full_tensor(@tensor, dim, poly_order);
num_basis = size(indx_mat, 1);
sample = importdata('CCsmolyak_d3lvl4_x.dat');
weight = importdata('CCsmolyak_d3lvl4_w.dat');
num_sample = length(weight);
measure_mat = zeros(num_sample, num_basis);

for k = 1:num_sample
  measure_mat(k,:) = eval_tensor_poly(@legendre_norm, sample(k,:)', dim, ...
                                      poly_order, indx_mat)';
end

training = f(sample);

coef = measure_mat'*diag(weight)*training;
