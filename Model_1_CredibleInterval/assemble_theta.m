function theta_full = assemble_theta(theta, para_fixed, para_change_exclude_idx)
% ASSEMBLE_THETA - Reconstruct full parameter vector by inserting fixed values
% 
% Usage:
%   theta_full = assemble_theta(theta, para_fixed, para_change_exclude_idx)
%
% Inputs:
%   theta                  - vector of variable (unknown) parameters
%   para_fixed             - vector of fixed parameter values
%   para_change_exclude_idx - indices (1-based) where para_fixed values should be inserted
%
% Output:
%   theta_full             - full parameter vector with both fixed and variable values

    % Total number of parameters
    total_num_params = length(theta) + length(para_fixed);

    % Initialize the full parameter vector
    theta_full = zeros(1, total_num_params);

    % Insert fixed parameters into specified indices
    theta_full(para_change_exclude_idx) = para_fixed;

    % Find remaining indices for variable parameters
    all_indices = 1:total_num_params;
    change_idx = setdiff(all_indices, para_change_exclude_idx);

    % Insert variable parameters
    theta_full(change_idx) = theta;
end