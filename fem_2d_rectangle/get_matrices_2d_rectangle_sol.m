function [M_R, Srr_R, Sss_R, Dr_R, Ds_R] = get_matrices_2d_rectangle_sol(k)
% get_matrices_2d_rectangle    Matrices for Poisson using FEM in 2D
%   get_matrices_2d_rectangle(k) generates the mass matrix M_R, the 
%   stiffness matrices Srr_R and Sss_R, and the differentiation matrices 
%   Dr_R and Ds_R for continuous k-th order polynomial approximations on 
%   the reference rectangle.
%
%   Parameters:
%     - k : polynomial order for the approximate solution
%
%   Returns:
%     - M_R : Mass matrix on the reference interval
%     - Srr_R : Stiffness matrix on the reference interval
%     - Sss_R : Stiffness matrix on the reference interval
%     - Dr_R : Differentiation matrix with respect to r on the reference 
%              rectangle
%     - Ds_R : Differentiation matrix with respect to s on the reference 
%              rectangle
end