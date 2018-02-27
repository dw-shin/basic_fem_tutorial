% Exercise 2. Modify fem_for_poisson_1d.m to solve the Poisson problem with
% non-homogeneous Dirichlet boundary condition,
%
% -u''(x) = f(x) in \Omega \\
% u(x) = u_D(x) on \partial\Omega. \\
%
% After you modify fem\_for\_poisson\_1d.m, run this code in command window.
% >> fem_1d_exercise2()
%
%
%
%
% The following is the main code. If your code is right, the convergence
% rate is k. 
% You can also compare the solution obtained from fem_for_poisson_1d.m 
% with that from fem_for_poisson_1d_ex2_sol.m
%
% k = 1;
% iter = 6;
% a = 0;
% b = 1;
% M = 2.^(1+(1:iter));
% f = @(x) 25*pi^2*sin(5*pi*x) - 2;
% u_D = @(x) x.^2;
% Du = @(x) 5*pi*cos(5*pi*x) + 2*x;
% 
% h = 1./M;
% count = zeros(3,1);
% error = zeros(1,iter);
% for j=1:iter
%     [c4n, n4e, n4db, ind4e] = mesh_fem_1d(a, b, M(j), k);
%     [M_R, S_R, D_R] = get_matrices_1d(k);
%     u = fem_for_poisson_1d(c4n, n4e, n4db, ind4e, M_R, S_R, f, u_D);
%     error(j) = compute_error_fem_1d(c4n, ind4e, M_R, D_R, u, Du);
% end
% convergence_rate = (log(error(2:end)) - log(error(1:end-1))) ./ (log(h(2:end)) - log(h(1:end-1)));
% str1 = sprintf('Convergence rate for k = %d \n', k);
% disp(str1)
% disp(num2str(convergence_rate))