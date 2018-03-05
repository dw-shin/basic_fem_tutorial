function [u, A, b, fns] = fem_for_poisson_2d_triangle_ex2_sol(c4n,n4e,n4db, ...
    ind4e,M_R,Srr_R,Srs_R,Ssr_R,Sss_R,f,u_D)
%fem_for_poisson_2d_triangle_ex2_sol    FEM solver for Poisson problem in 2D with
%                               triangular elements
%   fem_for_poisson_2d_triangle_ex2_sol(c4n,n4e,n4db,ind4e,M_R,Srr_R,Srs_R,Ssr_R, 
%   Sss_R,f,u_D) solves the Poisson problem. In order to use this code, 
%   mesh information (c4n,n4e,n4db,ind4e), matrices (M_R,Srr_R,Srs_R,Ssr_R, 
%   Sss_R), the source f, and the boundary condition u_D. Then the results 
%   of this code are the numerical solution u, the global stiffness matrix 
%   A, the global load vector b and the freenodes.
%
%   Parameters:
%     - c4n : coordinates for nodes.
%     - n4e : nodes for elements.
%     - n4db : nodes for Dirichlet boundary.
%     - ind4e : indices for elements
%     - M_R : Mass matrix on the reference triangle
%     - Srr_R : Stiffness matrix on the reference triangle
%     - Srs_R : Stiffness matrix on the reference triangle
%     - Ssr_R : Stiffness matrix on the reference triangle
%     - Sss_R : Stiffness matrix on the reference triangle
%     - f : RHS in the Poisson problem
%     - u_D : Dirichlet boundary condition for the solution u
%
%   Returns:
%     - u : numerical solution
%     - A : Global stiffness matrix
%     - b : Global right-hand side
%     - fns : free nodes
end