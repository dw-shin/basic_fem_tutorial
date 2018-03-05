function error = compute_error_fem_2d_triangle(c4n,n4e,ind4e,M_R,Dr_R,Ds_R,u,ux,uy)
%compute_error_fem_2d_triangle    Semi H1 error (2D triangular element)
%   compute_error_fem_2d_triangle(c4n,n4e,ind4e,M_R,Dr_R,Ds_R,u,ux,uy) 
%   computes the semi H1 error between the exact solution and the FE solution.
%
%   Parameters:
%     - c4n : coordinates for nodes.
%     - n4e : nodes for elements.
%     - ind4e : indices for elements
%     - M_R : Mass matrix on the reference triangle
%     - Dr_R : Differentiation matrix with respect to r on the reference 
%              triangle
%     - Ds_R : Differentiation matrix with respect to s on the reference 
%              triangle
%     - u : numerical solution
%     - ux : Derivative of the exact solution with respect to x for the 
%            model problem
%     - uy : Derivative of the exact solution with respect to y for the 
%            model problem
%
%   Returns:
%     - error  Semi H1 error between the exact solution and the FE solution.

error = 0;
for j=1:size(n4e,1)
    xr = (c4n(n4e(j,1),1)-c4n(n4e(j,3),1))/2; 
    yr = (c4n(n4e(j,1),2)-c4n(n4e(j,3),2))/2;
    xs = (c4n(n4e(j,2),1)-c4n(n4e(j,3),1))/2; 
    ys = (c4n(n4e(j,2),2)-c4n(n4e(j,3),2))/2;
    J = xr*ys-xs*yr;
    rx=ys/J; ry=-xs/J; sx=-yr/J; sy=xr/J;
    
    Dx_u = (rx*Dr_R+sx*Ds_R)*u(ind4e(j,:));
    Dy_u = (ry*Dr_R+sy*Ds_R)*u(ind4e(j,:));
    Dex=ux(c4n(ind4e(j,:),:)) - Dx_u;
    Dey=uy(c4n(ind4e(j,:),:)) - Dy_u;
    error=error+J*(Dex'*M_R*Dex+Dey'*M_R*Dey);
end
error=sqrt(error);
end