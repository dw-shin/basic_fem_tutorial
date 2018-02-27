iter = 10;
a = 0;      % left-end point of the domain
b = 1;      % right-end point of the domain
k = 3;      % polynomial order for the approximate solution
M = 2.^(1+(1:iter));        % the number of elements
f = @(x) 25*pi^2*sin(5*pi*x);        % RHS in the Poisson problem
u_D = @(x) x*0;     % Dirichlet boundary condition for the solution u
Du = @(x) 5*pi*cos(5*pi*x);     % Derivative of the exact solution for the model problem

error = zeros(1,iter);
h = 1./M;
for j=1:iter
    [c4n, n4e, n4db, ind4e] = mesh_fem_1d(a, b, M(j), k);    
    [M_R, S_R, D_R] = get_matrices_1d(k);
    u = fem_for_poisson_1d(c4n, n4e, n4db, ind4e, M_R, S_R, f, u_D);
    error(j) = compute_error_fem_1d(c4n, ind4e, M_R, D_R, u, Du);
end
rateE = (log(error(2:end)) - log(error(1:end-1))) ./ (log(h(2:end)) - log(h(1:end-1)));
disp(rateE)
