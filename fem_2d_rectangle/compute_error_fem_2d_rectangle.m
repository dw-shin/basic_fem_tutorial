function error = compute_error_fem_2d_rectangle(c4n, n4e, ind4e, M_R, ...
    Dr_R, Ds_R, u, ux, uy)
% compute_error_fem_2d_rectangle    Semi H1 error (2D rectangular element)
%   compute_error_fem_2d_rectangle(c4n,n4e,ind4e,M_R,Dr_R,Ds_R,u,ux,uy) 
%   computes the semi H1 error between the exact solution and the FE solution.
%
%   Parameters:
%     - c4n : coordinates for nodes.
%     - n4e : nodes for elements.
%     - ind4e : indices for elements
%     - M_R : Mass matrix on the reference interval
%     - Dr_R : Differentiation matrix along r-direction on the reference 
%              interval
%     - Ds_R : Differentiation matrix along s-direction on the reference 
%              interval
%     - u : numerical solution
%     - ux : Derivative of the exact solution with respect to x for the 
%            model problem
%     - uy : Derivative of the exact solution with respect to y for the 
%            model problem
%
%   Returns:
%     - error : Semi H1 error between the exact solution and the FE solution.

error = 0;
for j=1:size(ind4e,1)
    xr = (c4n(n4e(j,2),1)-c4n(n4e(j,1),1))/2;
    ys = (c4n(n4e(j,4),2)-c4n(n4e(j,1),2))/2;
    J = xr*ys;
    rx = ys/J; sy = xr/J;
    
    Dx_u = rx*Dr_R*u(ind4e(j,:));
    Dy_u = sy*Ds_R*u(ind4e(j,:));
    Dex = ux(c4n(ind4e(j,:),:)) - Dx_u;
    Dey = uy(c4n(ind4e(j,:),:)) - Dy_u;
    error = error + J*(Dex'*M_R*Dex + Dey'*M_R*Dey);
end
error = sqrt(error);
end