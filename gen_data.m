% Generate input and output samples by using a function as the black box. In
% practice, the black box is your simulation code, experiment, etc.
%
% Author : Xiu Yang
% Date   : 12/20/2016

function [] = gen_data(dim, num_sample)
  % "Black box", here we use a third-order polynomial
  f = @(x) sum(x) + 0.25*sum(x)^2 + 0.025*sum(x)^3;
  % Generate input samples
  rand('seed', 1);
  % The input is a MxP matrix, where M is the number of samples and P is the
  % dimension.
  input = rand(dim, num_sample)'*2-1;
  % Compute the output samples from the black box
  output = zeros(num_sample, 1);
  for k = 1:num_sample
    output(k) = f(input(k,:));
  end

  col_pnt = importdata('uniform_d12_lvl3_x.dat');
  num_col_pnt = size(col_pnt, 1);
  col_output = zeros(num_col_pnt, 1);
  for k = 1:num_col_pnt
    col_output(k) = f(col_pnt(k,:));
  end
  
  save('data.mat', 'input', 'output', 'col_output');
