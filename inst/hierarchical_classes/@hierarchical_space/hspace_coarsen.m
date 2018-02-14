% HSPACE_COARSEN: coarsen the hierarchical space, updating the fields of the object.
%
%   [hspace, Ccoarse] = hspace_coarsen (hspace, hmsh, funs_to_reactivate, removed_cells, reactivated_cell)
%
% INPUT:
%
%   hspace:             object representing the fine hierarchical space (see hierarchical_space)
%   hmsh:               object representing the hierarchical mesh, already coarsened (see hierarchical_mesh)
%   funs_to_reactivate: cell array with the indices, in the tensor product setting, of the functions to reactivate for each level
%   removed_cells:      cell array with the elements removed during coarsening, for each level
%   reactivated_cell:   cell array with the elements to be reactivated
%
% OUTPUT:
%
%   hspace:      object representing the coarsened hierarchical space (see hierarchical_space)
%   u_coarse:    coarsened dofs
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017 Luca Coradello, Rafael Vazquez
%
% Copyright (C) 2017 Massimo Carraturo
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [hspace, u] = hspace_coarsen (hspace, hmsh, FtR, removed_cells, reactivated_cell)

boundary = ~isempty (hspace.boundary);

% Update active functions
if(nargout < 2 && ~isempty(find(~cellfun(@isempty,removed_cells),1)))
    hspace = update_active_functions (hspace, hmsh, FtR, reactivated_cell, removed_cells);
else
    [hspace, u_coarse] = update_active_functions (hspace, hmsh, FtR, reactivated_cell, removed_cells);
    if(isempty(find(~cellfun(@isempty,removed_cells),1)))
        u = hspace.dofs;
    end
end

% Update the matrices for changing basis
hspace.Csub = hspace_subdivision_matrix (hspace, hmsh);
        
% Reconstruct the hierachical dofs structure
ndof_per_level = cellfun (@numel, hspace.active);

if(nargout >= 2 && hspace.truncated && ~isempty(find(~cellfun(@isempty,removed_cells),1)))
    ndof_lev = cumsum (ndof_per_level(1:hspace.nlevels));
    u = zeros(hspace.ndof, 1);
    % loop over levels
    for lev = hspace.nlevels:-1:2
        % if the level has active dofs
        if ~isempty(hspace.active{lev})
            u(ndof_lev(lev-1)+1:ndof_lev(lev)) =  hspace.Csub{lev}(:,ndof_lev(lev-1)+1:ndof_lev(lev))'*u_coarse{lev};
        end
    end
    u(1:ndof_lev(1)) = hspace.Csub{1}(:,1:ndof_lev(1))'*u_coarse{1};
elseif (nargout >= 2 && ~hspace.truncated && ~isempty(find(~cellfun(@isempty,removed_cells),1)))
    u = zeros(sum (ndof_per_level), 1);
    weights = cell(1, hspace.nlevels);
    for lev = hspace.nlevels:-1:1
        weights{lev} = cell(1, max(hspace.active{lev}));
        for iFunc = 1:max(hspace.active{lev})
            weights{lev}{iFunc} = cell(1, hspace.nlevels);
            for iiLev = hspace.nlevels:-1:1
                weights{lev}{iFunc}{iiLev} = cell(1, max(hmsh.active{iiLev}));
            end
        end
    end
    
    % loop over levels from n to 1
    for lev = hspace.nlevels:-1:1
        ndof_until_lev = sum (ndof_per_level(1:lev-1));

        [~, funs_support] = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), hspace.active{lev});
        for iFun=1:numel(hspace.active{lev})
            % function supports activated in projection
            index2smoothVec = funs_support{iFun};
            hier_func_index = ndof_until_lev+find(ismember(hspace.active{lev},hspace.active{lev}(iFun)));

            % loop over function supports
            for j=1:numel(index2smoothVec)
                % if active
                if ~isempty(find(ismember(index2smoothVec(j),hmsh.active{lev}),1))
                    c_den = 0.0;
                    % loop over control values...
                    [C, hier_basis_functions] = mlbzrextr(hmsh, hspace, hspace.active, lev, index2smoothVec(j));
                    c_num = sum(C(hier_basis_functions==hier_func_index, :));
                    
                    % ... loop over each support element
                    for i_el=1:numel(index2smoothVec)
                        [C, hier_basis_functions] = mlbzrextr(hmsh, hspace, hspace.active, lev, index2smoothVec(i_el));
                        
                        c_den_i = sum(C(hier_basis_functions==hier_func_index, :));
                        c_den = c_den_i + c_den;
                    end
                    
                    weights{lev}{hspace.active{lev}(iFun)}{lev}{index2smoothVec(j)} = c_num/c_den;
                else
                    children_cells = hmsh_get_active_children(hmsh, lev, index2smoothVec(j));
                    
                    for u_level = 1:numel(children_cells)
                        if ~isempty(children_cells{u_level})
                            for iChild = 1:length(children_cells{u_level})
                                c_den = 0.0;
                                [C, hier_basis_functions] = mlbzrextr(hmsh, hspace, hspace.active, lev+u_level-1, children_cells{u_level}(iChild));
                                c_num = sum(C(hier_basis_functions==hier_func_index, :))/((2^hmsh.ndim)^(u_level-1));
                                for i_el=1:numel(index2smoothVec)
                                    if ~isempty(find(ismember(index2smoothVec(i_el),hmsh.active{lev}),1))
                                        [C, hier_basis_functions] = mlbzrextr(hmsh, hspace, hspace.active, lev, index2smoothVec(i_el));
                                        c_den_i = sum(C(hier_basis_functions==hier_func_index, :));
                                        c_den = c_den_i + c_den;
                                    else
                                        index2smooth_den = index2smoothVec(i_el);
                                        uchildren_cells = hmsh_get_active_children(hmsh, lev, index2smooth_den);
                                        for uu_level = 1:numel(uchildren_cells)
                                            for iiChild=1:length(uchildren_cells{uu_level})
                                                [C, hier_basis_functions] = mlbzrextr(hmsh, hspace, hspace.active, lev+uu_level-1, uchildren_cells{uu_level}(iiChild));
                                                c_den_i = sum(C(hier_basis_functions==hier_func_index, :))/((2^hmsh.ndim)^(uu_level-1));
                                                c_den = c_den_i + c_den;
                                            end
                                        end
                                    end
                                end
                                weights{lev}{hspace.active{lev}(iFun)}{lev+u_level-1}{children_cells{u_level}(iChild)} = c_num/c_den;
                            end
                        end
                    end                    
                end
            end
            % end function supports loop
            
        end
        for iDof = 1:length(hspace.active{lev})
            active_dof = hspace.active{lev}(iDof);
            
            hier_support = hspace_get_function_support(hspace, hmsh, hspace.active, lev, active_dof);
            for u_lev = 1:numel(hier_support)
                if ~isempty(hier_support{u_lev})
                    for iSupp = 1:length(hier_support{u_lev})
                        active_support = hier_support{u_lev}(iSupp);
                        [R, hierc_funct_indices] = mlrcnst(hmsh, hspace, hspace.active, hspace.deactivated, lev+u_lev-1, active_support);
                        hierc_values =  R'*u_coarse{lev+u_lev-1}{active_support};
                        
                        u(ndof_until_lev+find(ismember(hspace.active{lev},active_dof))) =...
                            weights{lev}{active_dof}{lev+u_lev-1}{active_support}*...
                            hierc_values(hierc_funct_indices==(ndof_until_lev+find(ismember(hspace.active{lev},active_dof))))+...
                            u(ndof_until_lev+find(ismember(hspace.active{lev},active_dof)));
                        
                    end
                end
            end
        end
    end
end

% Fill the information for the boundaries
if (boundary)
    Nf = cumsum ([0, hspace.ndof_per_level]);
    for iside = 1:2*hmsh.ndim
        if (hmsh.ndim > 1)
            FtR_boundary = cell (size (FtR));
            for lev = 1:numel (FtR)
                [~,~,FtR_boundary{lev}] = intersect (FtR{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
            end
            cells_boundary = cell (size (removed_cells));
            for lev = 1:numel (removed_cells)
                cells_boundary{lev} = get_boundary_indices (iside, hmsh.mesh_of_level(lev).nel_dir, removed_cells{lev});
            end
            reactivated_cells_boundary = cell (size (reactivated_cell));
            for lev = 1:numel (reactivated_cell)
                reactivated_cells_boundary{lev} = get_boundary_indices (iside, hmsh.mesh_of_level(lev).nel_dir, reactivated_cell{lev});
            end
            hspace.boundary(iside) = hspace_coarsen (hspace.boundary(iside), hmsh.boundary(iside),...
                FtR_boundary, cells_boundary, reactivated_cells_boundary);
            
            nlevels_aux = hspace.boundary(iside).nlevels;
        elseif (hmsh.ndim == 1)
            nlevels_aux = hspace.nlevels;
        end
        
        dofs = [];
        for lev = 1:nlevels_aux
            [~,iact] = intersect (hspace.active{lev}, hspace.space_of_level(lev).boundary(iside).dofs);
            dofs = union (dofs, Nf(lev) + iact);
        end
        hspace.boundary(iside).dofs = dofs;
    end
    
else
    hspace.boundary = [];
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [hspace, u_coarse] = update_active_functions (hspace, hmsh, funs_to_reactivate, reactivated_cell, removed_cells)
%
% This function updates the active (hspace.active) and deactivated (hspace.deactivated) degrees of freedom,
% reactivating the functions in marked_funs.
% The function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:                     the fine space, an object of the class hierarchical_space
%           hmsh:                       an object of the class hierarchical_mesh, already coarsened
%           funs_to_reactivate:         indices of active functions of level lev to be reactivated
%           reactivated_cell:           indices of cells of level lev to be reactivated
%           removed_cells:              cell array with the elements removed during coarsening, for each level
%
% Output:   hspace:    the coarsened space, an object of the class hierarchical_space
%           u_coarse:  coarsened_dofs
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
% Copyright (C) 2017 Massimo Carraturo
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [hspace, u_coarse] = update_active_functions (hspace, hmsh, funs_to_reactivate, reactivated_cell, removed_cells)

% active = hspace.active;
% deactivated = hspace.deactivated;
% active_funs_support = cell(1, hspace.nlevels);                      % cell array with the active cells supporting reactivated functions
% Bzr_ext_container = cell(hmsh.ndim, hspace.nlevels);                % n-dimensional cell array to cache Bézier extraction matrices of each level
%
% % initialize the cell array to store the coarsened dofs u_coarse
% u_coarse = cell(1, hspace.nlevels);
% ndof_per_level_prev = cellfun (@numel, hspace.active);
% ndof_until_lev = sum (ndof_per_level_prev(1:hspace.nlevels-1));
% ndof_above_lev = sum (ndof_per_level_prev(1:hspace.nlevels));
% u_coarse{hspace.nlevels} = hspace.Csub{hspace.nlevels}(:,ndof_until_lev+1:ndof_above_lev)*hspace.dofs(ndof_until_lev+1:ndof_above_lev);
%
% funs2remove = sp_get_basis_functions (hspace.space_of_level(hspace.nlevels), hmsh.mesh_of_level(hspace.nlevels), removed_cells{hspace.nlevels});
% active{hspace.nlevels} = setdiff(active{hspace.nlevels}, funs2remove);
%
% for lev = hspace.nlevels-1:-1:1
%     ndof_until_lev = sum (ndof_per_level_prev(1:lev-1));
%     ndof_above_lev = sum (ndof_per_level_prev(1:lev));
%
%     u_coarse{lev} = hspace.Csub{lev}(:,ndof_until_lev+1:ndof_above_lev)*hspace.dofs(ndof_until_lev+1:ndof_above_lev);
%
%     [reactivated_funs_support, ~] = sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), funs_to_reactivate{lev});
%
%     funs2remove = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), removed_cells{lev});
%     active{lev} = setdiff(active{lev}, funs2remove);
%     active{lev} = union (active{lev}, funs_to_reactivate{lev});
%
%     deactivated{lev} = setdiff (deactivated{lev}, funs_to_reactivate{lev});
%
%     % active cells in projection
%     active_reactivated_funs_support = intersect(hmsh.active{lev}, reactivated_funs_support);         % active cells of the support of reactivated functions
%     active_funs_support{lev} = union(reactivated_cell{lev}, active_reactivated_funs_support);        % cells active in projection
%
% end
active = hspace.active;
deactivated = hspace.deactivated;
active_funs_support = cell(1, hspace.nlevels);                      % cell array with the active cells supporting reactivated functions
Bzr_ext_container = cell(hmsh.ndim, hspace.nlevels);                % n-dimensional cell array to cache Bézier extraction matrices of each level
reactivated_funs_support = cell(1, hspace.nlevels);

u_coarse = cell(1, hspace.nlevels);
ndofs_per_level = cellfun (@numel, hspace.active);
ndof_until_lev = sum (ndofs_per_level(1:hspace.nlevels-1));
ndof_above_lev = sum (ndofs_per_level(1:hspace.nlevels));
u_coarse{hspace.nlevels} = hspace.Csub{hspace.nlevels}(:,ndof_until_lev+1:ndof_above_lev)*hspace.dofs(ndof_until_lev+1:ndof_above_lev);


for lev = hspace.nlevels:-1:2
    ndof_until_lev = sum (ndofs_per_level(1:lev-2));
    ndof_above_lev = sum (ndofs_per_level(1:lev-1));
    
    u_coarse{lev-1} = hspace.Csub{lev-1}(:,ndof_until_lev+1:ndof_above_lev)*hspace.dofs(ndof_until_lev+1:ndof_above_lev);
    
    reactivated_funs_support{lev-1} = sp_get_cells (hspace.space_of_level(lev-1), hmsh.mesh_of_level(lev-1), funs_to_reactivate{lev-1});
    
    if (strcmpi (hspace.type, 'standard') && ~isempty (removed_cells{lev}))
        removed_funs = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), removed_cells{lev});
        active{lev} = setdiff (active{lev}, removed_funs);
        active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
        deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
        
    elseif (strcmpi (hspace.type, 'simplified') && ~isempty (funs_to_reactivate{lev-1}))
        active{lev-1} = union (active{lev-1}, funs_to_reactivate{lev-1});
        deactivated{lev-1} = setdiff (deactivated{lev-1}, funs_to_reactivate{lev-1});
        children = hspace_get_children (hspace, lev-1, funs_to_reactivate{lev-1});
        
        neighbors = sp_get_neighbors (hspace.space_of_level(lev-1), hmsh.mesh_of_level(lev-1), funs_to_reactivate{lev-1});
        deact_neighs = intersect (deactivated{lev-1}, neighbors);
        children = setdiff (children, hspace_get_children (hspace, lev-1, deact_neighs));
        active{lev} = setdiff (active{lev}, children);
    end
    
    % active cells in projection
    active_reactivated_funs_support = intersect(hmsh.active{lev-1}, reactivated_funs_support{lev-1});         % active cells of the support of reactivated functions
    active_funs_support{lev-1} = union(reactivated_cell{lev-1}, active_reactivated_funs_support);      % cells active in projection
    
end

if (nargout >= 2 && ~isempty(find(~cellfun(@isempty,removed_cells),1)))
    if (~hspace.truncated)
        clear u_coarse
        u_coarse_temp = cell(1, hspace.nlevels-1);
        
        % loop over levels from n to 1
        for lev = hspace.nlevels:-1:1

            u_coarse_temp{lev} = cell(1, max(hmsh.active{lev}));
            cells_activated_in_proj = hmsh.active{lev}; %sp_get_cells (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), function_activated_in_proj);
            % loop over elements involved in projection
            
            for el = 1:numel(cells_activated_in_proj)
                % if the element was already active
                if isempty(find(ismember(cells_activated_in_proj(el),reactivated_cell{lev}),1)) &&...
                        ~isempty(find(ismember(cells_activated_in_proj(el),hmsh.active{lev}),1))
                    [C, hier_funct] = mlbzrextr(hmsh, hspace, hspace.active, lev, cells_activated_in_proj(el));
                    u_coarse_temp{lev}{cells_activated_in_proj(el)} = C' * hspace.dofs(hier_funct);
                % if reactivated
                elseif ~isempty(find(ismember(cells_activated_in_proj(el),reactivated_cell{lev}),1))
                    % get children of the cell
                    children_cells = hmsh_get_children(hmsh, lev, cells_activated_in_proj(el));
                    [i,j,k] = ind2sub(hmsh.mesh_of_level(lev).nel_dir, cells_activated_in_proj(el));
                    bzr_values = 0;
                    [I,J,K] = ind2sub(hmsh.mesh_of_level(lev+1).nel_dir, children_cells);
                    % project children onto lev
                    for iChild=1:length(children_cells)
                        [C_child, hier_funct] = mlbzrextr(hmsh, hspace, hspace.active, lev+1, children_cells(iChild));
                        child_sub = [I(iChild),J(iChild),K(iChild)];
                        M = bzrproj_tp(hspace.space_of_level(lev+1).degree, hmsh, lev+1, lev, child_sub, [i,j,k]);
                        bzr_values = M * C_child' * hspace.dofs(hier_funct) + bzr_values;
                    end
                    % reconstruct hierarchical functions
                    u_coarse_temp{lev}{cells_activated_in_proj(el)} = bzr_values;
                % if deactivated but support of activated function in
                % projection
                else
                    [C, hier_funct] = mlbzrextr(hmsh, hspace, hspace.active, lev, hmsh.active{lev}(el));
                    u_coarse_temp{lev}{hmsh.active{lev}(el)} = C' * hspace.dofs(hier_funct);
                    
                end
            end
            % end loop over active functions of level
        end
        % end loop over levels      
        u_coarse = u_coarse_temp;        
        
        %% Computation of the matrix to pass from the original to the coarsened space
    else
        % Bèzier extraction of the finest level
        knots = hspace.space_of_level(hspace.nlevels).knots;
        degree = hspace.space_of_level(hspace.nlevels).degree;
        nel_dir = cell(1, hspace.nlevels);
        nel_dir{hspace.nlevels} = zeros(1, hmsh.ndim);
        % loop over dimensions
        for idim=1:hmsh.ndim
            Bzr_ext_container{idim, hspace.nlevels} = bzrextr(knots{idim}, degree(idim));
            nel_dir{hspace.nlevels}(idim) = numel(knots{idim})-2*degree(idim)-1;
        end
        
        % loop over levels from n-1 to 1
        for lev = hspace.nlevels-1:-1:1
            % Bèzier extraction of level lev
            knots = hspace.space_of_level(lev).knots;
            for idim=1:hmsh.ndim
                Bzr_ext_container{idim, lev} = bzrextr(knots{idim}, degree(idim));
                nel_dir{lev}(idim) = numel(knots{idim})-2*degree(idim)-1;
            end
            % if level contains elements to be coarsened...
            if  ~isempty (reactivated_cell{lev})
                % express active_funs_supported as linear combination of funs of
                % level lev+1
                ndof_above_lev = sum (ndofs_per_level(1:lev+1));
                % express lower levels funs as linear combination of funs of
                % level lev
                active_dofs_lc = hspace.Csub{lev+1}(:,1:ndof_above_lev)*hspace.dofs(1:ndof_above_lev);
                % get cells to be activated in projection
                cells_activated_in_proj = active_funs_support{lev};
                % initialize temporary coarse dofs vector
                u_coarse_temp = cell(1, max(cells_activated_in_proj));
                % element dimensional scaling parameter
                htarget_el = hmsh.msh_lev{lev}.element_size(1);
                
                % loop over elements involved in projection
                for el = 1:numel(cells_activated_in_proj)
                    % get children of the cell
                    children_cells = hmsh_get_children(hmsh, lev, cells_activated_in_proj(el));
                    % get indeces of functions to be activated in projection
                    funs_projected = sp_get_basis_functions (hspace.space_of_level(lev), hmsh.mesh_of_level(lev), cells_activated_in_proj(el));
                    % initialize u_coarse_temp of the element el
                    u_coarse_temp{cells_activated_in_proj(el)} = zeros(hspace.space_of_level(lev).ndof, 1);
                    % functions with support on children cells to be
                    % projected
                    funs_supported = sp_get_basis_functions (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), children_cells);
                    % initialize projection tensor
                    C = 1;
                    B_extr_temp = 1;
                    % get indeces for projection
                    [I,J,K] = ind2sub(nel_dir{lev+1}, children_cells);
                    children_cells_unidim_indeces = [I, J, K];
                    % check size of cells_activated_in_proj
                    if (size(cells_activated_in_proj,2) < size(cells_activated_in_proj, 1))
                        cells_activated_in_proj = cells_activated_in_proj';
                    end
                    [I,J,K] = ind2sub(nel_dir{lev}, cells_activated_in_proj);
                    cells_activated_in_proj_dir = [I; J; K];
                    
                    % loop over dimensions
                    for idim = 1:hmsh.ndim
                        hsource_el = htarget_el/hmsh.nsub(idim);
                        source_cell_indeces = unique(children_cells_unidim_indeces(:,idim));
                        % get local Bézier projector
                        B_el_proj = bzrproj_el_hcoarse( hspace.space_of_level(lev).degree(idim),...
                            Bzr_ext_container{idim, lev+1}, inv(Bzr_ext_container{idim, lev}(:,:,cells_activated_in_proj_dir(idim,el))),...
                            hsource_el, htarget_el, source_cell_indeces);
                        % assembly the operators of each sub-element
                        B_target_el = zeros(hspace.space_of_level(lev).degree(idim)+1, numel(source_cell_indeces) + hspace.space_of_level(lev).degree(idim));
                        for j=1:numel(source_cell_indeces)
                            B_target_el(:,j:j+hspace.space_of_level(lev).degree(idim)) = B_target_el(:,j:j+hspace.space_of_level(lev).degree(idim)) + B_el_proj(:,:,j);
                        end
                        % kronecker product
                        C = kron (B_target_el, C);
                        B_extr_temp = kron (Bzr_ext_container{idim, lev}(:,:,cells_activated_in_proj_dir(idim,el)), B_extr_temp);
                    end
                    % end loop over dimensions
                    Bzr_ext(:,:,cells_activated_in_proj(el)) = B_extr_temp;
                    
                    %project
                    u_coarse_temp{cells_activated_in_proj(el)}(funs_projected) = u_coarse_temp{cells_activated_in_proj(el)}(funs_projected) + C*active_dofs_lc(funs_supported);
                end
                % end elements loop
                
                % smoothing of active functions with support on reactiveted cells
                u_coarse{lev}(funs_to_reactivate{lev}) = smooth_dofs (hspace, hmsh, u_coarse_temp, funs_to_reactivate{lev}, lev, Bzr_ext);
            end
            % end if
        end
        % end loop over levels
    end
end

hspace.active = active(1:hspace.nlevels);
hspace.deactivated = deactivated(1:hspace.nlevels);
hspace.ndof_per_level = cellfun (@numel, hspace.active);
hspace.ndof = sum (hspace.ndof_per_level);

hspace.coeff_pou = ones (hspace.ndof, 1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function smoothed_dofs = smooth_dofs (hspace, hmsh, u_coarse, funs_to_smooth, level)
%
% This function updates the active (hspace.active) and deactivated (hspace.deactivated) degrees of freedom,
% reactivating the functions in marked_funs.
% The function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
%
% Input:    hspace:                                the fine space, an object of the class hierarchical_space
%           hmsh:                                  an object of the class hierarchical_mesh, already coarsened
%           u_coarse [1 x nel]:                    cell array of projected dofs over each target element
%           funs_to_smooth:                        functions to be averaged
%           level:                                 hierarchical mesh level
%
% Output:   smoothed_dofs:              array of smoothed dofs
%
% Copyright (C) 2017 Massimo Carraturo
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
function smoothed_dofs = smooth_dofs (hspace, hmsh, u_coarse, funs_to_smooth, level, C)

% functions supports
[~, funs_support_index] = sp_get_cells (hspace.space_of_level(level), hmsh.mesh_of_level(level), funs_to_smooth);

% initialize smoothed dofs vector
smoothed_dofs = zeros(max(funs_to_smooth), 1);
% loop over functions to smooth
for i=1:numel(funs_to_smooth)
    % cell supports activated in projection
    index2smooth{1} = intersect(funs_support_index{i}, find(~cellfun('isempty', u_coarse)));
    % indeces of functions
    nel_funs_support = numel(index2smooth{1});
    fun_index = zeros(1, nel_funs_support);
    for j=1:nel_funs_support
        fun_index(j) = find(sp_get_basis_functions (hspace.space_of_level(level), hmsh.mesh_of_level(level), index2smooth{1}(j))==funs_to_smooth(i));
    end
    % weights for smoothing
    %     w = wgheval_LLSQ( index2smooth );
    % loop over cell supports
    for j=1:numel(index2smooth{1})
        w = wgheval( C, fun_index, index2smooth{1}(j), index2smooth{1});
        smoothed_dofs(funs_to_smooth(i)) = smoothed_dofs(funs_to_smooth(i)) + u_coarse{index2smooth{1}(j)}(funs_to_smooth(i))*w;
    end
    % end j loop
end
% end loop over funs to smooth
smoothed_dofs = smoothed_dofs(funs_to_smooth);
end
