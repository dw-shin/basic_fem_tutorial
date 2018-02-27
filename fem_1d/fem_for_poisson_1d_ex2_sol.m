function [u, A, b, freenodes] = fem_for_poisson_1d_ex2_sol(c4n, n4e, ...
    n4db, ind4e, M_R, S_R, f, u_D)
%fem_for_poisson_1d_ex2_sol    FEM solver for Poisson problem in 1D
%   fem_for_poisson_1d_ex3_sol(c4n,n4e,n4db,ind4e,M_R,S_R,f,u_D) solves the
%   Poisson problem. In order to use this code, mesh information (c4n, n4e,
%   n4db, ind4e), matrices (M_R, S_R), the source f, and the boundary 
%   condition u_D. Then the results of this code are the numerical solution
%   u, the global stiffness matrix A, the global load vector b and the 
%   freenodes.
%
%   Parameters:
%     - c4n : coordinates for All nodes.
%     - n4e : nodes for elements.
%     - n4db : nodes for Dirichlet boundary.
%     - ind4e : index for elements
%     - M_R : Mass matrix on the reference interval
%     - S_R : Stiffness matrix on the reference interval
%     - f : RHS in the Poisson problem
%     - u_D : Dirichlet boundary condition for the solution u
%
%   Returns:
%     - u : numerical solution
%     - A : Global stiffness matrix
%     - b : Global right-hand side
%     - freenodes : free nodes
end