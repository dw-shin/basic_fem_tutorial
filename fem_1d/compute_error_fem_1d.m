function error = compute_error_fem_1d(c4n, ind4e, M_R, D_R, u, Du)
% compute_error_fem_1d    Semi H1 error 
%   compute_error_fem_1d(c4n,ind4e,M_R,D_R,u,Du) computes the semi H1 error 
%   between the exact solution and the FE solution.
%
%   Parameters:
%     - c4n : coordinates for All nodes.
%     - n4db : nodes for Dirichlet boundary.
%     - M_R : Mass matrix on the reference interval
%     - D_R : Differentiation matrix on the reference interval
%     - u : numerical solution
%     - Du : Derivative of the exact solution for the model problem
%
%   Returns:
%     - error : Semi H1 error between the exact solution and the FE solution.

error = 0;
for j = 1:size(ind4e,1)
    J = (c4n(ind4e(j,end)) - c4n(ind4e(j,1))) / 2;
    u_x = D_R*u(ind4e(j,:)) / J;
    De = Du(c4n(ind4e(j,:))) - u_x;
    error = error + J*De'*M_R*De;
end
error = sqrt(error);
end