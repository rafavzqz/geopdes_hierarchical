function [A, W, R] = prepare_to_refine_space(space)
%
% function [A, W, R] = prepare_to_refine_space(space)
%

Nf = cumsum([0 space.ndof_per_level]); 

A = cell(space.nlevels+1,1); 
W = cell(space.nlevels+1,1);

% El siguiente loop seguramente se puede evitar usando mat2cell
for lev = 1:space.nlevels
    ind_f = (Nf(lev)+1):Nf(lev+1);
    A{lev} = space.globnum_active(ind_f, 2:end);
    W{lev} = space.coeff(ind_f);
    W{lev} = W{lev}(:);
end

R = space.removed;

A{space.nlevels+1} = zeros(0,space.ndim);
W{space.nlevels+1} = zeros(0,1);
R{space.nlevels+1} = zeros(0,space.ndim);