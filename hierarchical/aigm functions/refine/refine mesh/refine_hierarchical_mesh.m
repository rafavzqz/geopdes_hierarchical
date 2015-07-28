function [hmsh, new_elements] = refine_hierarchical_mesh(hmsh, M, indices, boundary)
%
% function [hmsh, new_elements] = refine_hierarchical_mesh(hmsh, M, indices, boundary)
%
% This function computes the information for updating a hierarchical mesh
% when enlarging the underlying subdomains with some marked elements.
% ATENCION: Luego voy a limpiar mas este codigo
%
% Input:        hmsh: struct for current the hierarchical mesh
%               M:  (cell array),
%                   where M{lev} is a matrix whose rows are the indices
%                   of marked elements of level lev, for lev = 1:hmsh.nlevels
%               indices: (cell array),
%                   where indices{lev} satisfies M{lev} =
%                   hmsh.active{lev}(indices{lev}).
%                   If indices is [], then this function computes this
%                   information.
%               boundary: true or false, default: true. (Fill the
%                   information for the boundaries of the mesh).
%
% Output:   hmsh: struct for the new hierarchical mesh after refinement
%           new_elements: (cell-array), new_elements{lev}
%               contains the global indices of the new active cells of level
%               lev after refinement.
%

% This function uses:     update_active_cells
%

% ATENCION: Modificar nsub en Computation of mesh_of_level(nlevels+1) para
% refinamiento no diadico


if nargin == 3
    boundary = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of mesh_of_level(nlevels+1) if a new level will be active
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(M{hmsh.nlevels}) % if a new level is activated
    %tic
    %disp('Enlarging mesh_of_level:')
    % We will fill hmsh.mesh_of_level(nlevels+1)
    rule = msh_gauss_nodes (hmsh.mesh_of_level(1).nqn_dir);
    nsub = 2 * ones (1, hmsh.ndim); 
    [aaa,zeta] = kntrefine (hmsh.mesh_of_level(hmsh.nlevels).breaks, nsub-1, ones(1, hmsh.ndim), zeros (1, hmsh.ndim));
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    hmsh.mesh_of_level(hmsh.nlevels+1) = msh_cartesian (zeta, qn, qw, hmsh.geometry,'boundary', boundary);
    %tempo = toc;
    %fprintf(' %f seconds\n', tempo);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of active elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%tic
%disp('Updating cells:')
[hmsh, new_cells] = update_active_cells(hmsh, M, indices);

if (boundary && hmsh.ndim > 1)
    new_elements.interior = new_cells;
else
    new_elements = new_cells;
end

%tempo = toc;
%fprintf(' %f seconds\n', tempo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of msh_lev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
disp('Computing msh_lev:')
hmsh.msh_lev = cell(hmsh.nlevels,1);
% Lo siguiente se puede mejorar para reducir el tiempo de calculo
for ilev = 1 : hmsh.nlevels
    hmsh.msh_lev{ilev} = msh_evaluate_element_list(hmsh.mesh_of_level(ilev), hmsh.active{ilev});
end
tempo = toc;
fprintf(' %f seconds\n', tempo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the information for the boundaries
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
            M_sub = cell(1,hmsh.ndim);
            [M_sub{:}] = ind2sub(hmsh.mesh_of_level(lev).nel_dir,  M{lev}(:));
            M_sub = cell2mat(M_sub);
            M_boundary{lev} = M_sub(M_sub(:,i) == boundary_ind(lev),ind);
            % Mejorar lo siguiente
            switch hmsh.boundary(iside).ndim
                case 2,
                    M_boundary{lev} = sub2ind(hmsh.boundary(iside).mesh_of_level(lev).nel_dir,M_boundary{lev}(:,1),M_boundary{lev}(:,2));
                case 3,
                    M_boundary{lev} = sub2ind(hmsh.boundary(iside).mesh_of_level(lev).nel_dir,M_boundary{lev}(:,1),M_boundary{lev}(:,2),M_boundary{lev}(:,3));
            end
            %end
        end
        [hmsh.boundary(iside), new_elements.boundary{iside}] = refine_hierarchical_mesh(hmsh.boundary(iside), M_boundary, [], false);
    end
else
    hmsh.boundary = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%