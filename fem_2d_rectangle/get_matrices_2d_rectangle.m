function [M_R, Srr_R, Sss_R, Dr_R, Ds_R] = get_matrices_2d_rectangle(k)
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

if k==1
    M_R = [4 2 2 1; 2 4 1 2; 2 1 4 2; 1 2 2 4]/9;
    Srr_R = [2 -2 1 -1; -2 2 -1 1; 1 -1 2 -2; -1 1 -2 2]/6;
    Sss_R = [2 1 -2 -1; 1 2 -1 -2; -2 -1 2 1; -1 -2 1 2]/6;
    Dr_R = [-1 1 0 0; -1 1 0 0; 0 0 -1 1; 0 0 -1 1]/2;
    Ds_R = [-1 0 1 0; 0 -1 0 1; -1 0 1 0; 0 -1 0 1]/2;
elseif k==2
    M_R = [16 8 -4 8 4 -2 -4 -2 1; 8 64 8 4 32 4 -2 -16 -2;
        -4 8 16 -2 4 8 1 -2 -4; 8 4 -2 64 32 -16 8 4 -2;
        4 32 4 32 256 32 4 32 4; -2 4 8 -16 32 64 -2 4 8;
        -4 -2 1 8 4 -2 16 8 -4; -2 -16 -2 4 32 4 8 64 8;
        1 -2 -4 -2 4 8 -4 8 16]/225;
    Srr_R = [28 -32 4 14 -16 2 -7 8 -1; -32 64 -32 -16 32 -16 8 -16 8;
        4 -32 28 2 -16 14 -1 8 -7; 14 -16 2 112 -128 16 14 -16 2;
        -16 32 -16 -128 256 -128 -16 32 -16; 2 -16 14 16 -128 112 2 -16 14;
        -7 8 -1 14 -16 2 28 -32 4; 8 -16 8 -16 32 -16 -32 64 -32;
        -1 8 -7 2 -16 14 4 -32 28]/90;
    Sss_R = [28 14 -7 -32 -16 8 4 2 -1; 14 112 14 -16 -128 -16 2 16 2;
        -7 14 28 8 -16 -32 -1 2 4; -32 -16 8 64 32 -16 -32 -16 8;
        -16 -128 -16 32 256 32 -16 -128 -16; 8 -16 -32 -16 32 64 8 -16 -32;
        4 2 -1 -32 -16 8 28 14 -7; 2 16 2 -16 -128 -16 14 112 14;
        -1 2 4 8 -16 -32 -7 14 28]/90;
    Dr_R = [-3 4 -1 0 0 0 0 0 0; -1 0 1 0 0 0 0 0 0; 1 -4 3 0 0 0 0 0 0;
        0 0 0 -3 4 -1 0 0 0; 0 0 0 -1 0 1 0 0 0; 0 0 0 1 -4 3 0 0 0;
        0 0 0 0 0 0 -3 4 -1; 0 0 0 0 0 0 -1 0 1; 0 0 0 0 0 0 1 -4 3]/2;
    Ds_R = [-3 0 0 4 0 0 -1 0 0; 0 -3 0 0 4 0 0 -1 0; 0 0 -3 0 0 4 0 0 -1;
        -1 0 0 0 0 0 1 0 0; 0 -1 0 0 0 0 0 1 0; 0 0 -1 0 0 0 0 0 1;
        1 0 0 -4 0 0 3 0 0; 0 1 0 0 -4 0 0 3 0; 0 0 1 0 0 -4 0 0 3]/2;
else
    M_R = 0; Srr_R = 0; Sss_R = 0; Dr_R = 0; Ds_R = 0;
end
end