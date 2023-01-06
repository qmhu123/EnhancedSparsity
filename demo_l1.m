% Legendre polynomials.
% Combination of rotation and L1.
% Ridge function test demonstration.
% 
% Author : mengqi hu
% Date   : 12/09/2019
% Update : 12/30/2022

clear;
close all;

addpath Tensor
addpath spgl1-1.9
addpath optimazations

% formatOut = 'yyyymmdd';
% today = datestr(now,formatOut);
% filename = ['Log_',today,'.txt'];
% diary(filename)

dim = 12;
num_all_sample = 20000;

% Generate data for the demonstration, in practice comment this line.

load('data.mat');
all_sample = input;
all_obser  = output;

% Collocation points for examining the error
col_pnt = importdata('uniform_d12_lvl3_x.dat');
col_wgt = importdata('uniform_d12_lvl3_w.dat');

load('l1pm.mat')

disp('finish reading data');

poly_order = 3;
num_basis = nchoosek(dim+poly_order, poly_order);
indx_mat = full_tensor(@tensor, dim, poly_order);
%load(['kernel_dim' num2str(dim) '_P' num2str(poly_order) '.mat']);

full_meas_mat = zeros(num_all_sample, num_basis);

for k = 1:num_all_sample
    full_meas_mat(k,:) = eval_tensor_poly(@legendre_norm, all_sample(k,:)', dim, ...
        poly_order, indx_mat)';
end

num_col_pnt = length(col_wgt);
col_meas_mat = zeros(num_col_pnt, num_basis);
rot_col_meas_mat = zeros(num_col_pnt, num_basis);

for k = 1:num_col_pnt
    col_meas_mat(k,:) = eval_tensor_poly(@legendre_norm, col_pnt(k,:)', dim, ...
        poly_order, indx_mat)';
end
f_norm = sqrt(dot(col_output.^2, col_wgt));

% Number of rotations
num_rot_iter = 3;

% Number of trials (if multiple trials are needed)
% num_trial = 100;

num_trial = 1;

tstart = tic;







corr = zeros(num_trial, num_rot_iter+1);


for num_sample = 100
% for num_sample = 100:20:180
                pm.lambda  = l1lambda((num_sample/20)-4);
                pm.delta  = l1delta((num_sample/20)-4);
    % Relative mean square error
    rmse = zeros(num_trial, num_rot_iter+1);
    % Coefficients
    recon_coef = zeros(num_basis, num_rot_iter+1, num_trial);
    
    sample_pnt = zeros(num_sample, dim);
    training = zeros(num_sample, 1);
    meas_mat = zeros(num_sample, num_basis);
    
    
    for t = 1:num_trial
%         disp([num_sample,t])
        %                 shift = (t-1)*num_sample;
        shift = (t*5)*num_sample+20;
        sample_pnt = all_sample(shift+1:shift+num_sample,:);
        meas_mat = full_meas_mat(shift+1:shift+num_sample,:);
        training = all_obser(shift+1:shift+num_sample);
        
        
        corr(t,1) = coherence(meas_mat);
        c = CS_L1_uncon_ADMM(meas_mat, training, pm);
        
        recon_coef(:,1,t) = c;
        %display('reweighted l1 validation error');
        numeric = col_meas_mat*c;
        rmse(t,1) = sqrt(dot((numeric-col_output).^2, col_wgt))/f_norm;
        
        % Rotate
        rotate_mat = eye(dim);
        rotate_pnt = sample_pnt;
        rot_col_pnt = col_pnt;
        
        for  iter = 1:num_rot_iter
            A = zeros(dim, num_sample);
            for k = 1:num_sample
                A(:,k) = eval_tensor_dpoly(@legendre_norm, @dlegendre_norm, rotate_pnt(k,:)', dim, poly_order, indx_mat)'*c;
            end
            
            [U,S,V]=svd(rotate_mat*A);
            
            rotate_mat = U;
            rotate_pnt = sample_pnt*U;
            rot_col_pnt = col_pnt*U;
            
            
            for k = 1:num_sample
                meas_mat(k,:) = eval_tensor_poly(@legendre_norm, rotate_pnt(k,:)', dim, ...
                    poly_order, indx_mat)';
            end            
            
            corr(t,iter+1) = coherence(meas_mat);
            c = CS_L1_uncon_ADMM(meas_mat, training, pm);
            
            recon_coef(:,iter+1,t) = c;
            %display('rotation error')
            for k = 1:num_col_pnt
                rot_col_meas_mat(k,:) = eval_tensor_poly(@legendre_norm, rot_col_pnt(k,:)', dim, ...
                    poly_order, indx_mat)';
            end
            numeric = rot_col_meas_mat*c;
            rmse(t,1+iter) = sqrt(dot((numeric-col_output).^2, col_wgt))/f_norm;
        end
    end
    toc(tstart)
end

tellapsed = toc(tstart);
disp(tellapsed)

beep