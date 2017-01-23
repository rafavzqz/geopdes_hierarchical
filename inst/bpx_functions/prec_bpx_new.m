function y = prec_bpx_new (x, bpx, n_level, iopt, dim)

int_dofs = bpx(n_level).int_dofs;

xx = zeros (bpx(n_level).ndof,1);
yy = xx;
xx(int_dofs) = x;

if (iopt == 1) % Jacobi solver
  for ii = 2:n_level
    xx_1 = bpx(ii).Qi'*xx;

    A = bpx(ii).Ai(bpx(ii).new_dofs,bpx(ii).new_dofs);
    yy_1 = spdiags (spdiags(A,0), 0, size(A,1), size(A,2)) \ xx_1;
    yy = yy + bpx(ii).Qi*yy_1;
  end
  
elseif (iopt == 2) % symmetric Gauss-Seidel solver
  for ii = 2:n_level
    xx_1 = bpx(ii).Qi'*xx;

    A = bpx(ii).Ai(bpx(ii).new_dofs,bpx(ii).new_dofs);
    yy_1 = sparse(tril (A)) * (spdiags(spdiags(A,0),0,size(A,1),size(A,2)) \ sparse(triu(A))) \xx_1;
    yy = yy + bpx(ii).Qi*yy_1;
  end
  
elseif (iopt == 4) % Richardson
  for ii = 2:n_level
    xx_1 = bpx(ii).Qi'*xx;

    h = 2 ^ (-ii);
    coeff = h ^ (2 - dim);

    yy_1 = coeff * xx_1;
    yy = yy + bpx(ii).Qi*yy_1;
  end
  
elseif (iopt == 5) % Exact solver everywhere, and Jacobi at the finest
  for ii = 2:n_level
    xx_1 = bpx(ii).Qi'*xx;

    A = bpx(ii).Ai(bpx(ii).new_dofs,bpx(ii).new_dofs);
    yy_1 = A \ xx_1;
    yy = yy + bpx(ii).Qi*yy_1;
  end
%   for ii = n_level
%     xx_1 = bpx(ii).Qi'*xx;
% 
%     A = bpx(ii).Ai(bpx(ii).new_dofs,bpx(ii).new_dofs);
%     yy_1 = spdiags (spdiags(A,0), 0, size(A,1), size(A,2)) \ xx_1;
%     yy = yy + bpx(ii).Qi*yy_1;
%   end
end
clear yy_1;

% Apply direct solver in the coarsest level
xx_0 = bpx(1).Qi'*xx;
yy_0 = bpx(1).Ai(bpx(1).new_dofs,bpx(1).new_dofs)\xx_0;
yy = yy + bpx(1).Qi*yy_0;

y = yy(int_dofs);

end