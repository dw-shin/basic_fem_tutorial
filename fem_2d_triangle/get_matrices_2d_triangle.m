function [M_R, Srr_R, Srs_R, Ssr_R, Sss_R, Dr_R, Ds_R] = get_matrices_2d_triangle(k)
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

if k==1
    M_R = [2 1 1; 1 2 1; 1 1 2]/6;
    Srr_R = [1 -1 0; -1 1 0; 0 0 0]/2;
    Srs_R = [1 0 -1; -1 0 1; 0 0 0]/2;
    Ssr_R = [1 -1 0; 0 0 0; -1 1 0]/2;
    Sss_R = [1 0 -1; 0 0 0; -1 0 1]/2;
    Dr_R = [-1 1 0; -1 1 0; -1 1 0]/2;
    Ds_R = [-1 0 1; -1 0 1; -1 0 1]/2;
elseif k==2
    M_R = [6 0 -1 0 -4 -1; 0 32 0 16 16 -4; -1 0 6 -4 0 -1;
        0 16 -4 32 16 0; -4 16 0 16 32 0; -1 -4 -1 0 0 6]/90;
    Srr_R = [3 -4 1 0 0 0; -4 8 -4 0 0 0; 1 -4 3 0 0 0;
        0 0 0 8 -8 0; 0 0 0 -8 8 0; 0 0 0 0 0 0]/6;
    Srs_R = [3 0 0 -4 0 1; -4 4 0 4 -4 0; 1 -4 0 0 4 -1;
        0 4 0 4 -4 -4; 0 -4 0 -4 4 4; 0 0 0 0 0 0]/6;
    Ssr_R = [3 -4 1 0 0 0; 0 4 -4 4 -4 0; 0 0 0 0 0 0;
        -4 4 0 4 -4 0; 0 -4 4 -4 4 0; 1 0 -1 -4 4 0]/6;
    Sss_R = [3 0 0 -4 0 1; 0 8 0 0 -8 0; 0 0 0 0 0 0;
        -4 0 0 8 0 -4; 0 -8 0 0 8 0; 1 0 0 -4 0 3]/6;
    Dr_R = [-3 4 -1 0 0 0; -1 0 1 0 0 0; 1 -4 3 0 0 0;
        -1 2 -1 -2 2 0; 1 -2 1 -2 2 0; 1 0 -1 -4 4 0]/2;
    Ds_R = [-3 0 0 4 0 -1; -1 -2 0 2 2 -1; 1 -4 0 0 4 -1;
        -1 0 0 0 0 1; 1 -2 0 -2 2 1; 1 0 0 -4 0 3]/2;
else
    M_R = 0; Srr_R = 0; Srs_R = 0; Ssr_R = 0; Sss_R = 0; Dr_R = 0; Ds_R = 0;
end
end