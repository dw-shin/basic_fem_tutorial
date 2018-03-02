function [u, A, b, freenodes] = fem_for_poisson_2d_rectangle_ex3_sol(c4n, n4e, ...
    n4db, n4nb, ind4e, M_R, Srr_R, Sss_R, M_R1D, f, u_D, u_N)
% fem_for_poisson_2d_rectangle_ex3_sol    FEM solver for Poisson problem in 2D with
%   rectangular elements
%
%   fem_for_poisson_2d_rectangle_ex3_sol(c4n,n4e,n4db,n4nb,ind4e,M_R, ...
%   Srr_R,Sss_R,M_R1D,f,u_D,u_N) solves the Poisson problem. In order to
%   use this code, mesh information (c4n, n4e, n4db, ind4e), matrices (M_R,
%   Srr_R, Sss_R), the source f, and the boundary condition u_D. Then the
%   results of this code are the numerical solution u, the global stiffness 
%   matrix A, the global load vector b and the freenodes.
%
%   Parameters:
%     - c4n : coordinates for nodes.
%     - n4e : nodes for elements.
%     - n4db : nodes for Dirichlet boundary.
%     - n4nb : nodes for Neumann boundary.
%     - ind4e : indices for elements
%     - M_R : Mass matrix on the reference rectangle
%     - Srr_R : Stiffness matrix on the reference rectangle
%     - Sss_R : Stiffness matrix on the reference rectangle
%     - M_R1D : Mass matrix on the reference interval in 1D
%     - f : RHS in the Poisson problem
%     - u_D : Dirichlet boundary condition for the solution u
%     - u_N : Neumann boundary condition for the solution u
%
%   Returns:
%     - u : numerical solution
%     - A : Global stiffness matrix
%     - b : Global right-hand side
%     - freenodes : free nodes
end