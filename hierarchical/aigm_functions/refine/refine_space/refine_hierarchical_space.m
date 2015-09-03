function hspace = refine_hierarchical_space (hmsh, hspace, M, flag, new_cells)
%
% function hspace = refine_hierarchical_space(hmsh, hspace, M, flag, new_cells, boundary)
%
% This function updates the information in hspace. The variable hmsh has
% already been updated by refine_hierarchical_mesh
% ATENCION: Luego voy a limpiar mas esta funcion
%
%
% Input:        hmsh:
%               hspace:
%               M:  (cell array),
%                   where M{lev} contains the indices
%                   of marked functions or elements of level lev, for lev =
%                   1,2,...
%               flag: 'functions' or 'elements'
%               new_cells: see refine_hierarchical_mesh
%               boundary: true or false, default: true. (Fill the
%                   information for the boundaries of the space).
%
%
% Output:   hspace updated
%

% This function uses:    compute_functions_to_deactivate
%                        update_active_functions
%                        compute_matrices_for_changing_basis
%

boundary = ~isempty (hspace.boundary);

% Computation of indices of functions of level lev that will become
% nonactive when removing the functions or elements in M{lev}
M = compute_functions_to_deactivate(hmsh, hspace, M, flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Enlarging space_of_level and Proj, if it is needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (numel(hspace.space_of_level) < hmsh.nlevels)
  msh_level = hmsh.mesh_of_level(hmsh.nlevels);
  degree = hspace.space_of_level(hmsh.nlevels-1).degree;
  knots = kntrefine (hspace.space_of_level(hmsh.nlevels-1).knots, hmsh.nsub-1, degree, degree-1);
  hspace.space_of_level(hmsh.nlevels) = sp_bspline (knots, degree, msh_level);
  coarse_space = hspace.space_of_level(hmsh.nlevels-1).constructor (msh_level);
end

if (size(hspace.Proj,1) < (hmsh.nlevels - 1))
  for idim = 1:hmsh.ndim
    degree = coarse_space.sp_univ(idim).degree;
    knt_coarse = coarse_space.sp_univ(idim).knots;
    knt_fine   = hspace.space_of_level(hmsh.nlevels).sp_univ(idim).knots;
    hspace.Proj{hmsh.nlevels-1,idim} = basiskntins (degree, knt_coarse, knt_fine);        
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of active functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
disp('Updating active dofs:')

if (boundary && hmsh.ndim > 1)
    new_elem = new_cells.interior;
else
    new_elem = new_cells;
end

% Deactivation of functions which have to be removed and activation of the
% new ones
hspace = update_active_functions(hspace, hmsh, new_elem, M);

tempo = toc;
fprintf(' %f seconds\n', tempo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Updating C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
%disp('Computing C:')
aux_active = hspace.active;
dif = hmsh.nlevels - hspace.nlevels;
if dif
    aux_active = vertcat(hspace.active,cell(dif,1));
end
ndof = zeros(1,hmsh.nlevels);
for l = 1:hmsh.nlevels
    ndof(l) = hspace.space_of_level(l).ndof;
end
hspace.C = compute_matrices_for_changing_basis(hmsh.nlevels, aux_active, ndof, hspace.Proj);
%tempo = toc;
%fprintf(' %f seconds\n', tempo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of sp_lev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing sp_lev:')
hspace.sp_lev = cell(hmsh.nlevels,1);
for ilev = 1 : hmsh.nlevels
    hspace.sp_lev{ilev} = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'gradient', true,'hessian', true);
end
tempo = toc;
fprintf(' %f seconds\n', tempo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the information for the boundaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (boundary && hmsh.ndim > 1)
    for iside = 1:2*hmsh.ndim
        %%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
        ind = setdiff (1:hmsh.ndim, ceil(iside/2));
        i = setdiff (1:hmsh.ndim, ind);
        if mod(iside,2) == 1
            boundary_ind = ones(hspace.nlevels,1);
        else
            boundary_ind = zeros(hspace.nlevels,1);
            for lev = 1:hspace.nlevels
                boundary_ind(lev) = hspace.space_of_level(lev).ndof_dir(i);
            end
        end
        M_boundary = cell(size(M));
        for lev = 1:numel(M)
            % if ~isempty(M{lev})
            M_sub = cell(1,hmsh.ndim);
            [M_sub{:}] = ind2sub(hspace.space_of_level(lev).ndof_dir,  M{lev}(:));
            M_sub = cell2mat(M_sub);
            M_boundary{lev} = M_sub(M_sub(:,i) == boundary_ind(lev),ind);
            % Mejorar lo siguiente
            switch hmsh.boundary(iside).ndim %XXXXXX CHECK IF THIS IS CORRECT
                case 2,
                    M_boundary{lev} = sub2ind(hspace.boundary(iside).space_of_level(lev).ndof_dir,M_boundary{lev}(:,1),M_boundary{lev}(:,2));
                case 3,
                    M_boundary{lev} = sub2ind(hspace.boundary(iside).space_of_level(lev).ndof_dir,M_boundary{lev}(:,1),M_boundary{lev}(:,2),M_boundary{lev}(:,3));
            end
            %end
        end
        hspace.boundary(iside) = refine_hierarchical_space(hmsh.boundary(iside), hspace.boundary(iside), ...
            M_boundary, 'functions', new_cells.boundary{iside});
        % Now, we fill hspace.boundary(iside).dofs
        globnum_active_boundary = [hspace.boundary(iside).globnum_active(:,1:i) boundary_ind(hspace.boundary(iside).globnum_active(:,1)) ...
            hspace.boundary(iside).globnum_active(:,(i+1):end)];
        [unos, hspace.boundary(iside).dofs] = ismember(globnum_active_boundary,hspace.globnum_active,'rows');
        if any(unos~=1)
            disp('Warning: Error when computing hspace.boundary().dofs')
            pause
        end
    end
else
    hspace.boundary = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%