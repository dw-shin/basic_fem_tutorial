function [M_R, S_R, D_R] = get_matrices_1d(k)
%get_matrices_1d    Matrices for Poisson using FEM in 1D
%   get_matrices_1d(k) generates the mass matrix M_R, the stiffness
%   matrix S_R and the differentiation matrix D_R for continuous k-th 
%   order polynomial approximations on the reference interval.
%
%   Parameters:
%     - k : polynomial order for the approximate solution
%
%   Returns:
%     - M_R : Mass matrix on the reference interval
%     - S_R : Stiffness matrix on the reference interval
%     - D_R : Differentiation matrix on the reference interval

if k==1
    M_R = [2 1; 1 2]/3;
    S_R = [1 -1; -1 1]/2;
    D_R = [-1 1; -1 1]/2;
elseif k==2
    M_R = [4 2 -1; 2 16 2; -1 2 4]/15;
    S_R = [7 -8 1; -8 16 -8; 1 -8 7]/6;
    D_R = [-3 4 -1; -1 0 1; 1 -4 3]/2;
else
    M_R = 0; S_R = 0; D_R = 0;
end
end