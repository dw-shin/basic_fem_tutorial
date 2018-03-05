function [M_R, Srr_R, Srs_R, Ssr_R, Sss_R, Dr_R, Ds_R] = get_matrices_2d_triangle_sol(k)
%get_matrices_2d_triangle    Matrices for Poisson using FEM in 2D
%   get_matrices_2d_triangle(k) generates the mass matrix M_R, the 
%   stiffness matrices Srr_R, Srs_R, Ssr_R and Sss_R, and the 
%   differentiation matrices Dr_R and Ds_R for continuous k-th order 
%   polynomial approximations on the reference interval.
%
%   Parameters:
%     - k : polynomial order for the approximate solution
%
%   Returns:
%     - M_R : Mass matrix on the reference interval
%     - Srr_R : Stiffness matrix on the reference interval
%     - Srs_R : Stiffness matrix on the reference interval
%     - Ssr_R : Stiffness matrix on the reference interval
%     - Sss_R : Stiffness matrix on the reference interval
%     - Dr_R : Differentiation matrix with respect to r on the reference 
%              triangle
%     - Ds_R : Differentiation matrix with respect to s on the reference 
%              triangle
end