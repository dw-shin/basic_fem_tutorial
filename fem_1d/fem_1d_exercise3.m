% Exercise 3. Modify fem_for_poisson_1d_ex3.m to solve the Poisson problem with
% mixed boundary condition,
%
% -u''(x) = f(x) in \Omega \\
% u(x) = u_D(x) on \Gamma_D \\
% u'(x)n = u_N(x) on \Gamma_N.
%
% where \Gamma_D denotes the Dirichlet boundary, \Gamma_N denotes the 
% Neumann boundary, and n is the outward unit normal vector.
%
% After you modify fem_for_poisson_1d_ex3.m, run this code in command window.
% >> fem_1d_exercise3()
%
%
%
% The following is the main code. If your code is right, the convergence
% rate is k. 
% You can also compare the solution obtained from fem_for_poisson_1d_ex3.m 
% with that from fem_for_poisson_1d_ex3_sol.m
%
% k = 1;
% iter = 6;
% a = 0;
% b = 1;
% M = 2.^(1+(1:iter));
% f = @(x) 25*pi^2*sin(5*pi*x) - 2;
% u_D = @(x) x.^2;
% u_N = @(x) 5*pi*cos(5*pi*x) + 2*x;
% Du = @(x) 5*pi*cos(5*pi*x) + 2*x;
% 
% h = 1./M;
% error = zeros(1,iter);
% for j=1:iter
%     [c4n, n4e, n4db, ind4e] = mesh_fem_1d(a, b, M(j), k);
%     n4nb = n4db(2);
%     n4db = n4db(1);
%     [M_R, S_R, D_R] = get_matrices_1d(k);
%     u = fem_for_poisson_1d_ex3(c4n, n4e, n4db, n4nb, ind4e, M_R, S_R, f, u_D, u_N);
%     error(j) = compute_error_fem_1d(c4n, ind4e, M_R, D_R, u, Du);
% end
% convergence_rate = (log(error(2:end)) - log(error(1:end-1))) ./ (log(h(2:end)) - log(h(1:end-1)));
% str1 = sprintf('Convergence rate for k = %d \n', k);
% disp(str1)
% disp(num2str(convergence_rate))