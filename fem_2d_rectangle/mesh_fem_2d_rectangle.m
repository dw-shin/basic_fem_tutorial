function [c4n, n4e, ind4e, inddb] = mesh_fem_2d_rectangle(xl, xr, yl, yr, Mx, My, k)
% mesh_fem_2d_rectangle    Mesh geometry on 2D rectangular domain
%   mesh_fem_2d_rectangle(xl, xr, yl, yr, Mx, My, k) generates an 
%   uniform mesh on the domain [xl,xr]x[yl,yr] in 2D with Mx elements
%   along x-direction and My elements along y-direction. Also this code 
%   returns an index matrix for continuous k-th order polynomial 
%   approximations.
%
%   Parameters:
%     - xl : x-coordinate of bottom-left vertex of the domain
%     - xr : x-coordinate of top-right vertex of the domain
%     - yl : y-coordinate of bottom-left vertex of the domain
%     - yr : y-coordinate of top-right vertex of the domain
%     - Mx : the number of elements along x-direction
%     - My : the number of elements along y-direction
%     - k : polynomial order for the approximate solution
%
%   Returns:
%     - c4n : coordinates for nodes.
%     - n4e : nodes for elements.
%     - n4db : nodes for Dirichlet boundary.
%     - ind4e : indices for elements

% index
ind4e = zeros(Mx*My, (k+1)^2);
tmp = repmat((1:k:k*Mx)', 1, My) ...
    + repmat((0:k*(k*Mx+1):((k*Mx+1)*((My-1)*k+1)-1)), Mx, 1);
tmp = tmp(:);
for j=1:k+1
    ind4e(:,(j-1)*(k+1)+(1:(k+1))) = repmat(tmp+(j-1)*(k*Mx+1), 1, k+1) ...
        + repmat((0:k), Mx*My, 1);
end

% n4e
n4e = ind4e(:,[1 k+1 (k+1)^2 (k*(k+1)+1)]);

% indDb
inddb = [1:(k*Mx+1), 2*(k*Mx+1):(k*Mx+1):(k*Mx+1)*(k*My+1), ...
                ((k*Mx+1)*(k*My+1)-1):-1:(k*My*(k*Mx+1)+1), ...
                ((k*My-1)*(k*Mx+1)+1):-(k*Mx+1):(k*Mx+2)]';

% c4n
x = linspace(xl, xr, k*Mx+1);
y = linspace(yl, yr, k*My+1);
y = repmat(y, k*Mx+1, 1);
x = repmat(x, k*My+1, 1)';
c4n = [x(:), y(:)];
end