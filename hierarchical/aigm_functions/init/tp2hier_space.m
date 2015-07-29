function hspace = tp2hier_space (hmsh, space, space_type, boundary)
%
% function hspace = tp2hier_space (hmsh, space, space_type, boundary)
%
% This function initializes a struct hspace from the tensor product space
% and the hierarchical mesh hmsh already initializated with tp2hier_msh
%
% INPUT
%                       hmsh: see tp2hier_msh
%                       space: the spline space for the coarsest space (level 1)
%                       space_type: 0 (simplified basis), 1 (full basis)
%                       boundary: true or false. (Fill the information for
%                       the boundaries of the mesh).
%
% OUPUT
%                       struct hmsh (Hierarchical mesh)
% 
% FIELD_NAME    TYPE               DESCRIPTION
% ndim           (scalar)           number of parametric directions
% degree        (1 x ndim array)    polynomial degree in each parametric direction
% ncomp
% type          0 (simplified basis), 1 (full basis)
% ndof          (scalar)           total number of active functions 
% nlevels       (scalar)           the number of levels
% space_of_level (1 x nlevels)                tensor product space of each level, with 1d evaluations on the mesh of the same level (see sp_bspline)
% Proj           (hmsh.nlevels-1 x ndim cell-array) the coefficients relating 1D splines of two consecutive levels
%     Proj{l,i} is a matrix of size N_{l+1} x N_l where N_l is the number of univariate functions of level l in the direction l, such that
% A function B_{k,l} = \sum_j c^k_j B_{j,l+1}, and c_j = Proj{l,i}(j,k)
% globnum_active (ndof x (ndim+1))  global tensor-product numbering of active functions and their corresponding level
% ndof_per_level (1 x nlevels array) number of active functions on each level
% active        (1 x nlevels cell-array) List of active functions on each level
% coeff         (ndof x 1)         coefficientes to form the partition of the unity in the hierarchical space
% deactivated       (1 x nlevels cell-array) List of deactivated functions on each level
% C             (1 x hmsh.nlevels cell-array) Matrices for changing basis (see compute_matrices_for_changing_basis.m)
% sp_lev        (hmsh.nlevels x 1 cell-array) sp_lev{ilev} is a structure
% boundary
%

hspace.ndim = hmsh.ndim;
hspace.degree = space.degree;
hspace.ncomp = space.ncomp;
hspace.type = space_type;

hspace.nlevels = 1;
hspace.ndof = space.ndof;
hspace.active{1} = (1:space.ndof)';

aux = cell(hspace.ndim,1);
[aux{:}] = ind2sub(space.ndof_dir,1:space.ndof);
hspace.globnum_active = [ones(space.ndof,1) cell2mat(aux)'];

hspace.ndof_per_level = [space.ndof];
hspace.coeff = ones(space.ndof,1);    
hspace.deactivated{1} = zeros(0,1);  
hspace.space_of_level = space;
hspace.Proj = [];
hspace.C{1} = speye(space.ndof);
hspace.sp_lev{1} = sp_evaluate_element_list (hspace.space_of_level(1), hmsh.msh_lev{1}, 'gradient', true,'hessian', true);


if (boundary && hmsh.ndim > 1)
    for iside = 1:2*hmsh.ndim
        %%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
        ind = setdiff (1:hmsh.ndim, ceil(iside/2));
        i = setdiff (1:hmsh.ndim, ind);
        if mod(iside,2) == 1
            boundary_ind = 1;
        else
            boundary_ind = hspace.space_of_level(1).ndof_dir(i);
        end
        bound = tp2hier_space (hmsh.boundary(iside), space.boundary(iside), space_type, false);
        % Now, we fill hspace.boundary(iside).dofs
        globnum_active_boundary = [bound.globnum_active(:,1:i) boundary_ind*ones(bound.ndof,1) ...
            bound.globnum_active(:,(i+1):end)];    
        [unos, bound.dofs] = ismember(globnum_active_boundary,hspace.globnum_active,'rows');
        if any(unos~=1)
            disp('Warning: Error when computing hspace.boundary().dofs')
            pause
        end
        hspace.boundary(iside) = bound;
    end
else
    hspace.boundary = [];
end