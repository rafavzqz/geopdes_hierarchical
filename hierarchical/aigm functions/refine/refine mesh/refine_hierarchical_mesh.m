function [hmsh, new_elements] = refine_hierarchical_mesh(hmsh, M, indices, boundary)
%
% function [hmsh, new_elements] = refine_hierarchical_mesh(hmsh, M, indices, boundary)
%
% This function computes the information for updating a hierarchical mesh
% when enlarging the underlying subdomains with some marked elements.
% ATENCION: Luego voy a limpiar este codigo
%
% Input:        hmsh: struct for current the hierarchical mesh
%               M:  (1 x msh.nlevels cell array),
%                   where M{lev} is a matrix whose rows are the tensor-product indices
%                   of marked elements of level lev, for lev = 1:hmsh.nlevels
%               indices: (1 x msh.nlevels cell array),
%                   where indices{lev} satisfies M{lev} =
%                   hmsh.active{lev}(indices{lev}).
%                   If indices is [], then this function computes this
%                   information.
%               boundary: true or false, default: true. (Fill the
%                   information for the boundaries of the mesh).
%
% Output:   hmsh: struct for the new hierarchical mesh after refinement
%           new_elements: (msh.nlevels x 1 cell-array), new_elements{lev}
%               contains the global indices of the new active cells of level
%               lev after refinement.
%

% This function uses:     update_active_cells
%

if nargin == 3
    boundary = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of mesh_of_level(nlevels+1) if a new level will be active
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(M{hmsh.nlevels}) % if a new level is activated
    tic
    disp('Enlarging mesh_of_level:')
    nlevels = hmsh.nlevels + 1;
    nel_per_level = zeros(1, nlevels);
    % We will fill hmsh.mesh_of_level(nlevels+1)
    rule = msh_gauss_nodes (hmsh.mesh_of_level(1).nqn_dir);
    nsub = 2 * ones (1, hmsh.ndim);
    [aaa,zeta] = kntrefine (hmsh.mesh_of_level(hmsh.nlevels).breaks, nsub-1, ones(1, hmsh.ndim), zeros (1, hmsh.ndim));
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    hmsh.mesh_of_level(hmsh.nlevels+1) = msh_cartesian (zeta, qn, qw, hmsh.geometry,'boundary', boundary);
    tempo = toc;
    fprintf(' %f seconds\n', tempo);
else
    nlevels = hmsh.nlevels;
    nel_per_level = zeros(1, hmsh.nlevels);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We fill E with the information of active cells
Ne = cumsum([0; hmsh.nel_per_level(:)]);
E = cell(hmsh.nlevels,1);
% El siguiente loop se puede evitar usando mat2cell
for lev = 1:hmsh.nlevels
        ind_e = (Ne(lev)+1):Ne(lev+1);
        E{lev} = hmsh.globnum_active(ind_e, 2:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of active elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
disp('Updating cells:')
[E, hmsh.removed,new_cells] = update_active_cells(E, hmsh.removed, M, indices);

hmsh.nlevels = nlevels;
aux_ind_el = [];
for lev = 1:hmsh.nlevels
    nel_per_level(lev) = size(E{lev},1);
    aux_ind_el = vertcat(aux_ind_el, lev*ones(nel_per_level(lev),1));
end % for lev = 1:nlevels


if (boundary && hmsh.ndim > 1)
    new_elements.interior = new_cells;
else
    new_elements = new_cells;
end

hmsh.nel_per_level = nel_per_level;
hmsh.nel = sum(hmsh.nel_per_level);
hmsh.globnum_active = cell2mat(E);
hmsh.globnum_active = horzcat(aux_ind_el, hmsh.globnum_active);

hmsh.active = cell(hmsh.nlevels,1);
for ilev = 1:hmsh.nlevels
    % Mejorar lo siguiente
    switch hmsh.ndim
        case 1,
            hmsh.active{ilev} = E{ilev};
        case 2,
            hmsh.active{ilev} = sub2ind(hmsh.mesh_of_level(ilev).nel_dir,E{ilev}(:,1),E{ilev}(:,2));
        case 3,
            hmsh.active{ilev} = sub2ind(hmsh.mesh_of_level(ilev).nel_dir,E{ilev}(:,1),E{ilev}(:,2),E{ilev}(:,3));
    end
end
tempo = toc;
fprintf(' %f seconds\n', tempo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of msh_lev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
disp('Computing msh_lev:')
hmsh.msh_lev = cell(hmsh.nlevels,1);
% Unefficient way. This is ok for the first working version
for ilev = 1 : hmsh.nlevels
    elems = hmsh.active{ilev};
    hmsh.msh_lev{ilev} = msh_evaluate_element_list(hmsh.mesh_of_level(ilev), elems);
end
tempo = toc;
fprintf(' %f seconds\n', tempo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (boundary && hmsh.ndim > 1)
    for iside = 1:2*hmsh.ndim
        %%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
        ind = setdiff (1:hmsh.ndim, ceil(iside/2));
        i = setdiff (1:hmsh.ndim, ind);
        if mod(iside,2) == 1
            boundary_ind = ones(hmsh.nlevels,1);
        else
            boundary_ind = zeros(hmsh.nlevels,1);
            for lev = 1:hmsh.nlevels
                boundary_ind(lev) = hmsh.mesh_of_level(lev).nel_dir(i);
            end
        end
        M_boundary = cell(size(M));
        for lev = 1:numel(M)
           % if ~isempty(M{lev})
                M_boundary{lev} = M{lev}(M{lev}(:,i) == boundary_ind(lev),ind);
            %end
        end
        [hmsh.boundary(iside), new_elements.boundary{iside}] = refine_hierarchical_mesh(hmsh.boundary(iside), M_boundary, [], false);
    end
else
    hmsh.boundary = [];
end

