% Exercise 2. Modify fem_for_poisson_2d_rectangle_ex2.m to solve the Poisson 
% problem with non-homogeneous Dirichlet boundary condition,
%
% - div(grad u) = f(x) in \Omega \\
% u(x) = u_D(x) on \partial\Omega. \\
%
% After you modify fem_for_poisson_2d_rectangle_ex2.m, run this code in command window.
% >> fem_2d_rectangle_exercise2()
%
%
%
%
% The following is the main code. If your code is right, the convergence
% rate is k. 
% You can also compare the solution obtained from fem_for_poisson_2d_rectangle_ex2.m 
% with that from fem_for_poisson_2d_rectangle_ex2_sol.m
%
% iter = 5;
% k = 3;
% xl = 0; xr = 1; yl = 0; yr = 1; M = 2.^(1:iter);
% f=@(x) 2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)) - 2;
% u_D=@(x) x(:,1)*0 + x(:,1).^2;
% ux=@(x) pi*cos(pi*x(:,1)).*sin(pi*x(:,2)) + 2*x(:,1);
% uy=@(x) pi*sin(pi*x(:,1)).*cos(pi*x(:,2));
% 
% error=zeros(1,iter);
% time=zeros(1,iter);
% h=1./M;
% for j=1:iter
%     [c4n,n4e,ind4e,n4db] = mesh_fem_2d_rectangle(xl,xr,yl,yr,M(j),M(j),k);
%     [M_R, Srr_R, Sss_R, Dr_R, Ds_R] = get_matrices_2d_rectangle(k);
%     u = fem_for_poisson_2d_rectangle_ex2(c4n,n4e,n4db,ind4e,M_R,Srr_R,Sss_R,f,u_D);
%     error(j) = compute_error_fem_2d_rectangle(c4n,n4e,ind4e,M_R,Dr_R,Ds_R,u,ux,uy);
%     disp(j)
% end
% convergence_rate = (log(error(2:end)) - log(error(1:end-1))) ./ (log(h(2:end)) - log(h(1:end-1)));
% str1 = sprintf('Convergence rate for k = %d \n', k);
% disp(str1)
% disp(num2str(convergence_rate))