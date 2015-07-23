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
%               FIELD_NAME    TYPE               DESCRIPTION
%               ndim           (scalar)          number of parametric directions
%               rdim           (scalar)
%               nlevels       (scalar)           the number of levels
%               nsub          (1 x ndim array)           number of subdivisions on each level
%               mesh_of_level (1 x nlevels mesh) Cartesian mesh of each level (see msh_cartesian)
%               nel           (scalar)           total number of active cells  
%               nel_per_level (1 x nlevels array) number of active cells on each level
%               globnum_active (nel x (dim+1))  global tensor-product numbering of active cells and their corresponding level
%               active        (1 x nlevels cell-array) List of active elements on each level
%               removed       (1 x nlevels cell-array) List of removed cells on each level
%               msh_lev     (nlevels x 1 cell-array) msh_lev{ilev} is a structure
%               geometry
%               boundary    
%

hspace.ndim = hmsh.ndim;
hspace.degree = space.degree;
hspace.ncomp = space.ncomp;
hspace.type = space_type;

hspace.nlevels = 1;
hspace.ndof = space.ndof;
hspace.active{1} = 1:space.ndof;

aux = cell(hspace.ndim,1);
[aux{:}] = ind2sub(space.ndof_dir,1:space.ndof);
hspace.globnum_active = [ones(space.ndof,1) cell2mat(aux)'];

hspace.ndof_per_level = [space.ndof];
hspace.coeff = ones(space.ndof,1);    
hspace.removed{1} = zeros(0,hspace.ndim);  
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