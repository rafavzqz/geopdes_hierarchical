function [x, flag, relres, iter, resvec, eigest] = ...
                pcg_w_eigest( A, b, tol, maxit, M, x0, varargin )
%[X, FLAG, RELRES, ITER, RESVEC, EIGEST] = pcg( A, B, TOL, MAXIT, M, X0, ... )
%
%Solves the linear system of equations A*x = B by means of the  Preconditioned
%Conjugate Gradient iterative method.
%
%INPUT ARGUMENTS
%---------------
%
%A can be either a square (preferably sparse) matrix or a name of a function
%which computes A*x. In principle A should be symmetric and positive
%definite; if PCG finds A not positive definite, you will get a warning message
%and the FLAG output parameter will be set.
%
%B is the right hand side vector.
%
%TOL is the required relative tolerance for the residual error, B-A*X. The
%iteration stops if ||B-A*X|| <= TOL*||B-A*X0||, where ||.|| denotes the
%Euclidean norm. If TOL is empty or is omitted, the function sets TOL = 1e-6 by
%default.
%
%MAXIT is the maximum allowable number of iterations;  if [] is supplied for
%MAXIT, or PCG has less arguments, a default value equal to 20 is used.
%
%M is the (left) preconditioning matrix, so that the iteration is
%(theoretically) equivalent to solving by PCG P*x = M\B, with P = M\A.
%Note that a proper choice of the preconditioner may
%dramatically improve the overall performance of the method! Instead of matrix M, 
%the user may pass a function which returns the results of applying the inverse of
%M to a vector (usually this is the preferred way of using the preconditioner).
%If [] is supplied for M, or M is omitted, no preconditioning is
%applied.
%
%X0 is the initial guess. If X0 is empty or omitted, the function sets X0 to a
%zero vector by default.
%
%The arguments which follow X0 are treated as parameters, and passed in a proper
%way to any of the functions (A or M) which are passed to pcg. See the EXAMPLES
%for details.
%
%OUTPUT ARGUMENTS
%----------------
%
%X is the computed approximation to the solution x of A*x=B.
%
%FLAG reports on the convergence. FLAG = 0 means the solution converged
%and the tolerance criterion given by TOL is satisfied. FLAG = 1 means that the
%MAXIT limit for the iteration count was reached. FLAG = 3 reports that the
%(preconditioned) matrix was found not positive definite.
%
%RELRES is the ratio of the final residual to its initial value, measured in the
%Euclidean norm.
%
%ITER is the actual number of iterations performed.
%
%RESVEC describes the convergence history of the method. RESVEC(i,1) is the
%Euclidean norm of the residual, and RESVEC(i,2) is the preconditioned residual
%norm, after the (i-1)-th iteration, i = 1,2,...ITER+1. The preconditioned
%residual norm is defined as |||r|||^2 = r'*(M\r) where r = B-A*x, see also the
%description of M. If EIGEST is not required, only RESVEC(:,1) is returned.
%
%EIGEST returns the estimate for the smallest (EIGEST(1)) and largest
%(EIGEST(2)) eigenvalues of the preconditioned matrix P=M\A. 
%In particular, if no preconditioning is used, the
%estimates for the extreme eigenvalues of A are returned. EIGEST(1) is an
%overestimate and EIGEST(2) is an underestimate, so that EIGEST(2)/EIGEST(1) is
%a lower bound for cond(P,2), which nevertheless in the limit should
%theoretically be equal to the actual value of the condition number. 
%The method which computes EIGEST works only for symmetric positive
%definite A and M, and the user is responsible for verifying this assumption. 
%
%EXAMPLES 
%--------
%
%Let us consider a trivial problem with a diagonal matrix (we exploit the
%sparsity of A) 
%
%       N = 10; 
%       A = diag([1:N]); A = sparse(A);  
%       b = rand(N,1);
%
%EX. 1. Simplest use of PCG
%
%       x = pcg(A,b)
%
%EX. 2. PCG with a function which computes A*x
%
%       function y = applyA(x) y = [1:N]'.*x; end
%
%       x = pcg('applyA',b)
%
%EX. 3. Preconditioned iteration, with full diagnostics. The preconditioner (quite
%strange, because even the original matrix A is trivial) is defined as a
%function:
%
%       function y = applyM(x)          
%       K = floor(length(x)-2); 
%       y = x; 
%       y(1:K) = x(1:K)./[1:K]';        
%       end
%
%       [x, flag, relres, iter, resvec, eigest] = pcg(A,b,[],[],'applyM')
%       semilogy([1:iter+1], resvec);
%
%EX. 4. Finally, a preconditioner which depends on a parameter K:
%
%       function y = applyM(x, varargin)
%       K = varargin{1}; 
%       y = x; y(1:K) = x(1:K)./[1:K]';  
%       end
%
%       [x, flag, relres, iter, resvec, eigest] = pcg(A,b,[],[],'applyM',[],3)
%
%You can also run 
%
%       demo('pcg') 
%
%from the command line to see more simple examples of how the pcg works.
%
%SEE ALSO: sparse, pcr, gmres
%
%REFERENCES
%
%       [1] C.T.Kelley, 'Iterative methods for linear and nonlinear equations',
%       SIAM, 1995 (the base PCG algorithm) 
%       
%       [2] Y.Saad, 'Iterative methods for sparse linear systems', PWS 1996
%       (condition number estimate from PCG) Revised version of this book is
%       available online at http://www-users.cs.umn.edu/~saad/books.html
%

%% Copyright (C) 2004 Piotr Krzyzanowski <piotr.krzyzanowski@mimuw.edu.pl>
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%% 
%% REVISION HISTORY
%%
%% 2004-05-21, Piotr Krzyzanowski:
%%      Added 4 demos and 4 tests
%%
%% 2004-05-18, Piotr Krzyzanowski:
%%      Warnings use warning() function now
%%
%% 2004-04-29, Piotr Krzyzanowski:
%%      Added more warning messages when FLAG is not required
%%
%% 2004-04-28, Piotr Krzyzanowski:
%%      When eigest is required, resvec returns both the Euclidean and the
%%      preconditioned residual norm convergence history
%%
%% 2004-04-20, Piotr Krzyzanowski: 
%%      Corrected eigenvalue estimator. Changed the tridiagonal matrix
%%      eigenvalue solver to regular eig
%% 
if ((nargin < 6) || isempty(x0))
        x = zeros(size(b));
else
        x = x0;
end

if (nargin < 5)
        M = [];
end

if ((nargin < 4) || isempty(maxit))
        maxit = min(size(b,1),20);
end

maxit = maxit + 2;

if ((nargin < 3) || isempty(tol))
        tol = 1e-6;
end

preconditioned_residual_out = false;
if (nargout > 5)
        T = sparse(maxit,maxit);
        preconditioned_residual_out = true;
end

matrix_positive_definite = true;        % assume A is positive definite

if (~exist('OCTAVE_VERSION'))           % a trick to maintain source
                                        % compatibility with M*TLAB
        stderr = 2;
end

p = zeros(size(b)); 
oldtau = 1; 
if (isnumeric(A))                       % is A a matrix?
        r = b - A*x; 
else                                    % then A should be a function!
        r = b - feval(A,x,varargin{:});
end

%resvec(1,1) = norm(r); alpha = 1; iter = 2;

if (isnumeric(M))               % is M a matrix?
                if isempty(M)           %       if M is empty, use no precond
                        z = r;
                else                    %       otherwise, apply the precond
                        z = M\r;
                end
else                            % then M should be a function!
                z = feval(M,r,varargin{:});
end
resvec(1,1) = norm(z); alpha = 1; iter = 2;

while(  (resvec(iter-1,1) > tol*resvec(1,1)) && (iter < maxit) )
        
        if (isnumeric(M))               % is M a matrix?
                if isempty(M)           %       if M is empty, use no precond
                        z = r;
                else                    %       otherwise, apply the precond
                        z = M\r;
                end
        else                            % then M should be a function!
                z = feval(M,r,varargin{:});
        end
        tau = z'*r; resvec(iter-1,2) = sqrt(tau);
        beta = tau/oldtau;
        oldtau = tau;
        p = z + beta*p;
        if (isnumeric(A))               % is A a matrix?
                w = A*p;
        else                            % then A should be a function!
                w = feval(A,p,varargin{:});
        end
        oldalpha = alpha;                       % needed only for eigest
        alpha = tau/(p'*w);
        if (alpha <= 0.0) % negative matrix?
                matrix_positive_definite = false;
        end
        x = x + alpha*p;
        r = r - alpha*w;
        if (nargout > 5) && (iter>2)
                T(iter-1:iter, iter-1:iter) = T(iter-1:iter, iter-1:iter) + ...
                                [1 sqrt(beta); sqrt(beta) beta]./oldalpha;
                % EVS = eig(T(2:iter-1,2:iter-1));
                % fprintf(stderr,'PCG condest: %g (iteration: %d)\n', max(EVS)/min(EVS),iter);
        end;
%        resvec(iter,1) = norm(r);
        resvec(iter,1) = norm(z);
        iter = iter +1;
end


if (nargout > 5)
        if ( matrix_positive_definite )
                if (iter > 3)
                        T = T(2:iter-2,2:iter-2);
                        l = eig(T);
                        eigest = [min(l) max(l)];
                        % fprintf(stderr, 'PCG condest: %g\n',eigest(2)/eigest(1));
                else
                        eigest = [NaN NaN];
                        warning('PCG: eigenvalue estimate failed: iteration converged too fast.');
                end
        else
                eigest = [NaN NaN];
        end

        % apply the preconditioner once more and finish with the precond
        % residual
        if (isnumeric(M))               % is M a matrix?
                if isempty(M)           %       if M is empty, use no precond
                        z = r;
                else                    %       otherwise, apply the precond
                        z = M\r;
                end
        else                            % then M should be a function!
                z = feval(M,r,varargin{:});
        end
        resvec(iter-1,2) = sqrt(r'*z);
else
        resvec = resvec(:,1);   
end

flag = 0;
relres = resvec(iter-1,1)./resvec(1,1);
iter = iter - 2;
if ( iter >= (maxit-2) )
        flag = 1;
        if (nargout < 2)
                warning('PCG: maximum number of iterations (%d) reached\n',...
                        iter);
                warning('The initial residual norm was reduced %g times.\n',...
                        1.0/relres);
        end
else
        if (nargout < 2)
                fprintf(stderr, 'PCG: converged in %d iterations. ', iter);
                fprintf(stderr, 'The initial residual norm was reduced %g times.\n',...
                        1.0/relres);
        end     
end
if (~matrix_positive_definite)
        flag = 3;
        if (nargout < 2)
                warning('PCG: matrix not positive definite?\n');
        end
end

% end % durkbin babo
%!demo
%!
%!      # Simplest usage of pcg (see also 'help pcg')
%!
%!      N = 10; 
%!      A = diag([1:N]); b = rand(N,1); y = A\b; #y is the true solution
%!      x = pcg(A,b);
%!      printf('The solution relative error is %g\n', norm(x-y)/norm(y));
%!
%!      # You shouldn't be afraid if pcg issues some warning messages in this
%!      # example: watch out in the second example, why it takes N iterations 
%!      # of pcg to converge to (a very accurate, by the way) solution
%!demo
%!
%!      # Full output from pcg, except for the eigenvalue estimates
%!      # We use this output to plot the convergence history  
%!
%!      N = 10; 
%!      A = diag([1:N]); b = rand(N,1); X = A\b; #X is the true solution
%!      [x, flag, relres, iter, resvec] = pcg(A,b);
%!      printf('The solution relative error is %g\n', norm(x-X)/norm(X));
%!      title('Convergence history'); xlabel('Iteration'); ylabel('log(||b-Ax||/||b||)');
%!      semilogy([0:iter],resvec/resvec(1),'o-g;relative residual;');
%!demo
%!
%!      # Full output from pcg, including the eigenvalue estimates
%!      # Hilbert matrix is extremely ill conditioned, so pcg WILL have problems
%!
%!      N = 10; 
%!      A = hilb(N); b = rand(N,1); X = A\b; #X is the true solution
%!      [x, flag, relres, iter, resvec, eigest] = pcg(A,b,[],200);
%!      printf('The solution relative error is %g\n', norm(x-X)/norm(X));
%!      printf('Condition number estimate is %g\n', eigest(2)/eigest(1));
%!      printf('Actual condition number is   %g\n', cond(A));
%!      title('Convergence history'); xlabel('Iteration'); ylabel('log(||b-Ax||)');
%!      semilogy([0:iter],resvec,['o-g;absolute residual;';'+-r;absolute preconditioned residual;']);
%!demo
%!
%!      # Full output from pcg, including the eigenvalue estimates
%!      # We use the 1-D Laplacian matrix for A, and cond(A) = O(N^2)
%!      # and that's the reasone we need some preconditioner; here we take
%!      # a very simple and not powerful Jacobi preconditioner, 
%!      # which is the diagonal of A
%!
%!      N = 100; 
%!      A = zeros(N,N);
%!      for i=1:N-1 # form 1-D Laplacian matrix
%!              A(i:i+1,i:i+1) = [2 -1; -1 2];
%!      endfor
%!      b = rand(N,1); X = A\b; #X is the true solution
%!      maxit = 80;
%!      printf('System condition number is %g\n',cond(A));
%!      # No preconditioner: the convergence is very slow!
%!
%!      [x, flag, relres, iter, resvec, eigest] = pcg(A,b,[],maxit);
%!      printf('System condition number estimate is %g\n',eigest(2)/eigest(1));
%!      title('Convergence history'); xlabel('Iteration'); ylabel('log(||b-Ax||)');
%!      semilogy([0:iter],resvec(:,1),'o-g;NO preconditioning: absolute residual;');
%!
%!      pause(1);
%!      # Test Jacobi preconditioner: it will not help much!!!
%!
%!      M = diag(diag(A)); # Jacobi preconditioner
%!      [x, flag, relres, iter, resvec, eigest] = pcg(A,b,[],maxit,M);
%!      printf('JACOBI preconditioned system condition number estimate is %g\n',eigest(2)/eigest(1));
%!      hold on;
%!      semilogy([0:iter],resvec(:,1),'o-r;JACOBI preconditioner: absolute residual;');
%!
%!      pause(1);
%!      # Test nonoverlapping block Jacobi preconditioner: it will help much!
%!
%!      M = zeros(N,N);k=4
%!      for i=1:k:N # form 1-D Laplacian matrix
%!              M(i:i+k-1,i:i+k-1) = A(i:i+k-1,i:i+k-1);
%!      endfor
%!      [x, flag, relres, iter, resvec, eigest] = pcg(A,b,[],maxit,M);
%!      printf('BLOCK JACOBI preconditioned system condition number estimate is %g\n',eigest(2)/eigest(1));
%!      semilogy([0:iter],resvec(:,1),'o-b;BLOCK JACOBI preconditioner: absolute residual;');
%!      hold off;
%!test
%!
%!      #solve small diagonal system
%!
%!      N = 10; 
%!      A = diag([1:N]); b = rand(N,1); X = A\b; #X is the true solution
%!      [x, flag] = pcg(A,b,[],N+1);
%!      assert(norm(x-X)/norm(X),0,1e-10);
%!      assert(flag,0);
%!
%!test
%!
%!      #solve small indefinite diagonal system
%!      #despite A is indefinite, the iteration continues and converges
%!      #indefiniteness of A is detected
%!
%!      N = 10; 
%!      A = diag([1:N].*(-ones(1,N).^2)); b = rand(N,1); X = A\b; #X is the true solution
%!      [x, flag] = pcg(A,b,[],N+1);
%!      assert(norm(x-X)/norm(X),0,1e-10);
%!      assert(flag,3);
%!
%!test
%!
%!      #solve tridiagonal system, do not converge in default 20 iterations
%!
%!      N = 100; 
%!      A = zeros(N,N);
%!      for i=1:N-1 # form 1-D Laplacian matrix
%!              A(i:i+1,i:i+1) = [2 -1; -1 2];
%!      endfor
%!      b = ones(N,1); X = A\b; #X is the true solution
%!      [x, flag, relres, iter, resvec, eigest] = pcg(A,b,1e-12);
%!      assert(flag);
%!      assert(relres>1.0);
%!      assert(iter,20); #should perform max allowable default number of iterations
%!
%!test
%!
%!      #solve tridiagonal system with 'prefect' preconditioner
%!      #converges in one iteration, so the eigest does not work
%!      #and issues a warning
%!
%!      N = 100; 
%!      A = zeros(N,N);
%!      for i=1:N-1 # form 1-D Laplacian matrix
%!              A(i:i+1,i:i+1) = [2 -1; -1 2];
%!      endfor
%!      b = ones(N,1); X = A\b; #X is the true solution
%!      [x, flag, relres, iter, resvec, eigest] = pcg(A,b,[],[],A,b);
%!      assert(norm(x-X)/norm(X),0,1e-6);
%!      assert(flag,0);
%!      assert(iter,1); #should converge in one iteration
%!      assert(isnan(eigest),isnan([NaN NaN]));
%!
