function [c4n, n4e, n4db, ind4e] = mesh_fem_1d(a, b, M, k)
%mesh_fem_1d    Mesh geometry on 1 dimensional domain
%   mesh_fem_1d(a, b, M, k) generates an uniform mesh on the domain [a,b] 
%   in 1D with mesh size h = 1/M. Also this code returns an index matrix 
%   for continuous k-th order polynomial approximations.
%
%   Parameters:
%     - a : left-end point for the domain
%     - b : right-end point for the domain
%     - M : the number of elements
%     - k : polynomial order for the approximate solution
%
%   Returns
%     - c4n : coordinates for All nodes.
%     - n4e : nodes for elements.
%     - n4db : nodes for Dirichlet boundary.
%     - ind4e : index for elements
 
nrNodes = k*M + 1; % the number of nodes on the mesh in terms of k and N
c4n = linspace(a, b, nrNodes)'; 
n4e = [1:k:(nrNodes-1); (k+1):k:nrNodes]';
n4db = [1, nrNodes]';
ind4e = repmat(n4e(:,1), 1, k+1) + repmat((0:k), M, 1);
end