function [ hmsh ] = hmsh_balancing( hmsh, n_rb )
%HMSH_BALANCING relax the mesh avoiding two level jumps at the element
%interface
%
%     [ hmshb ] = hmsh_balancing( hmsh, n_rb )
%
% INPUT:
%
%     hmsh: the hierarchical mesh (see hierarchical_mesh)
%     n_rb: ring balancing paramenter
%
% OUTPUT:
%
%     hmshb: the balanced hierarchical mesh (see hierarchical_mesh)
%
% Copyright (C) 2017 Massimo Carraturo
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

old_elements = hmsh.active;
[hmsh, new_elements] = balance(hmsh, n_rb);
hmsh.msh_lev = update_msh_lev (hmsh, old_elements, new_elements);

% Update the boundary , calling the function recursively
boundary = ~isempty (hmsh.boundary);

if (boundary)
    if (hmsh.ndim > 1)
        for iside = 1:2*hmsh.ndim
            hmsh.boundary(iside) = balance_boundaries(hmsh.boundary(iside), n_rb);
        end
    end
else
    hmsh.boundary = [];
end


% Update hmsh.nel_per_level and hmsh.nel
hmsh.nel_per_level = cellfun (@numel, hmsh.active);
hmsh.nel = sum (hmsh.nel_per_level);

end

function [hmsh, new_elements] = balance(hmsh, n_rb)
new_elements = cell (hmsh.nlevels, 1);
active = hmsh.active;
deactivated = hmsh.deactivated;
%union of active and deactivated elements of hierarchical mesh
active_deactivated = cellfun(@union,active,deactivated,'UniformOutput',false);
el_dir = zeros(3, 1);
%loop over level
for lev=hmsh.nlevels:-1:2
    %loop over active elements of the level lev
    for el=1:numel(active{lev})
        %get parents of the active elements
        [p, ~] = hmsh_get_parent(hmsh, lev, active{lev}(el));
        %get neighbours of the parent
        % first get the local indeces in each direction ...
        [el_dir(1), el_dir(2), el_dir(3)] = ind2sub(hmsh.mesh_of_level(lev-1).nel_dir, p);
        % ... then get the parent's neighbors
        neighbors  = get_neighbors( hmsh, lev-1, el_dir(1), el_dir(2), el_dir(3), n_rb);
        inactive_neighbors = setdiff( neighbors, active_deactivated{lev-1});
        %set of parent of element e, p, and r_rb rings of neighbors of p
        %that lie around p
        P1 = union( neighbors, p);
        if (lev > 2)
            %get parents of all elements of p and cach them in P2
            P2 = [];
            for i=1:P1.size
                gp = hmsh_get_parent(hmsh, lev-1, P1(i));
                P2 = [P2 gp];
            end
            inactive_gp = setdiff( P2, active_deactivated{lev-2});
            %get all children of P2 and add them to P1
            children_P2 = [];
            for ii=1:P2.size
                children_P2 = [children_P2 hmsh_get_children(hmsh, lev-1, P2(ii))];
            end
            inactive_childrenP2 = setdiff( children_P2, active_deactivated{lev-1});
            %update the active elements in lev-1 and the deactivated elements
            %in l-2
            P1_inactive = union(inactive_neighbors, inactive_childrenP2);
            new_elements{lev-1} = union(new_elements{lev-1}, P1_inactive);
            hmsh.active{lev-1} = union(active{lev-1}, P1_inactive);
            hmsh.deactivated{lev-2} = union(deactivated{lev-2}, inactive_gp);
        else
            %update the active elements in lev-1
            hmsh.active{lev-1} = union(active{lev-1}, inactive_neighbors);
        end
    end
    
end
end

function  neighbors  = get_neighbors( hmsh, lev, ind_x, ind_y, ind_z, n_rb )
% x direction
neighbors_x_dir = linspace(ind_x-n_rb, ind_x+n_rb, 2*n_rb+1);
neighbors_x = zeros(size(neighbors_x_dir));
for i=1:numel(neighbors_x_dir)
    if (neighbors_x_dir(i)<=hmsh.mesh_of_level(lev).nel_dir(1) && neighbors_x_dir(i)>0)
        neighbors_x(i) = sub2ind(hmsh.mesh_of_level(lev).nel_dir, neighbors_x_dir(i), ind_y, ind_z);
    end
end
neighbors_x = unique(neighbors_x(neighbors_x>0));
% y direction
neighbors_y = [];
if numel(hmsh.mesh_of_level(lev).nel_dir) > 1
    neighbours_y_dir = linspace(ind_y-n_rb, ind_y+n_rb, 2*n_rb+1);
    neighbors_y = zeros(size(neighbours_y_dir));
    for i=1:numel(neighbours_y_dir)
        if(neighbours_y_dir(i)<=hmsh.mesh_of_level(lev).nel_dir(2)  && neighbours_y_dir(i)>0)
            neighbors_y(i) = sub2ind(hmsh.mesh_of_level(lev).nel_dir, ind_x, neighbours_y_dir(i), ind_z);
        end
    end
end
neighbors_y = unique(neighbors_y(neighbors_y>0));
% z direction
neighbors_z = [];
if numel(hmsh.mesh_of_level(lev).nel_dir) > 2
    neighbours_z_dir = linspace(ind_z-n_rb, ind_z+n_rb, 2*n_rb+1);
    neighbors_z = zeros(size(neighbours_z_dir));
    for i=1:numel(neighbours_z_dir)
        if (neighbours_z_dir(i)<=hmsh.mesh_of_level(lev).nel_dir(3)  && neighbours_z_dir(i)>0)
            neighbors_z(i) = sub2ind(hmsh.mesh_of_level(lev).nel_dir, ind_x, ind_y, neighbours_z_dir(i));% fill an auxiliary list of elements
        end
    end
end
neighbors_z = unique(neighbors_z(neighbors_z>0));

neighbors = unique([neighbors_x neighbors_y neighbors_z]);
end


function hmsh = balance_boundaries( hmsh, n_rb)
active = hmsh.active;
deactivated = hmsh.deactivated;
%union of active and deactivated elements of hierarchical mesh
active_deactivated = cellfun(@union,active,deactivated,'UniformOutput',false);
el_dir = zeros(3, 1);
%loop over level
for lev=hmsh.nlevels:-1:2
    %loop over active elements of the level lev
    for el=1:numel(active{lev})
        %get parents of the active elements
        [p, ~] = hmsh_get_parent(hmsh, lev, active{lev}(el));
        %get neighbours of the parent
        % first get the local indeces in each direction ...
        [el_dir(1), el_dir(2), el_dir(3)] = ind2sub(hmsh.mesh_of_level(lev-1).nel_dir, p);
        % ... then get the parent's neighbors
        neighbors  = get_neighbors_boundaries( hmsh, lev-1, el_dir(1), el_dir(2), el_dir(3), n_rb);
        inactive_neighbors = setdiff( neighbors, active_deactivated{lev-1});
        %set of parent of element e, p, and r_rb rings of neighbors of p
        %that lie around p
        P1 = union( neighbors, p);
        if (lev > 2)
            %get parents of all elements of p and cach them in P2
            P2 = [];
            for i=1:P1.size
                gp = hmsh_get_parent(hmsh, lev-1, P1(i));
                P2 = [P2 gp];
            end
            inactive_gp = setdiff( P2, active_deactivated{lev-2});
            %get all children of P2 and add them to P1
            children_P2 = [];
            for ii=1:P2.size
                children_P2 = [children_P2 hmsh_get_children(hmsh, lev-1, P2(ii))];
            end
            inactive_childrenP2 = setdiff( children_P2, active_deactivated{lev-1});
            %update the active elements in lev-1 and the deactivated elements
            %in l-2
            P1_inactive = union(inactive_neighbors, inactive_childrenP2);
            hmsh.active{lev-1} = union(active{lev-1}, P1_inactive);
            hmsh.deactivated{lev-2} = union(deactivated{lev-2}, inactive_gp);
        else
            %update the active elements in lev-1
            hmsh.active{lev-1} = union(active{lev-1}, inactive_neighbors);
        end
    end
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function msh_lev = update_msh_lev (hmsh, old_elements, new_elements)
%
% function msh_lev = update_msh_lev (hmsh, old_elements, new_elements)
%
% Update the information in msh_lev, computing only the elements that have been added to the mesh
%
% INPUT
%     hmsh: object representing the fine hierarchical mesh (see hierarchical_mesh)
%     old_elements{lev}: active elements in the previous coarse mesh
%     new_elements{lev}: elements that have been added after refinement
%
% OUTPUT
%     msh_lev: the structures with the msh information for each level
%

msh_lev = cell (hmsh.nlevels, 1);

for lev = 1:hmsh.nlevels
    if (lev > numel (old_elements) || numel (old_elements{lev}) == 0)
        msh_lev{lev} = msh_evaluate_element_list (hmsh.mesh_of_level(lev), hmsh.active{lev});
    else
        [~, iold, iold_act] = intersect (old_elements{lev}, hmsh.active{lev});
        msh_lev{lev}.ndim = hmsh.ndim;
        msh_lev{lev}.rdim = hmsh.rdim;
        msh_lev{lev}.nel = hmsh.nel_per_level(lev);
        msh_lev{lev}.elem_list = hmsh.active{lev}(:)';
        msh_lev{lev}.nel_dir = hmsh.mesh_of_level(lev).nel_dir;
        msh_lev{lev}.nqn_dir = hmsh.mesh_of_level(lev).nqn_dir;
        msh_lev{lev}.nqn = hmsh.mesh_of_level(lev).nqn;
        
        if (isempty (new_elements{lev}))
            indices = iold_act;
            msh_new = struct ('quad_weights', [], 'geo_map', [], 'geo_map_jac', [], 'geo_map_der2', [], 'jacdet', [], 'element_size', []);
        else
            msh_new = msh_evaluate_element_list (hmsh.mesh_of_level(lev), new_elements{lev});
            [~, ~, inew_act] = intersect (new_elements{lev}, hmsh.active{lev});
            indices = [iold_act(:); inew_act(:)];
        end
        msh_lev{lev}.quad_weights(:,indices) = [hmsh.msh_lev{lev}.quad_weights(:,iold), msh_new.quad_weights];
        msh_lev{lev}.geo_map(:,:,indices) = cat (3, hmsh.msh_lev{lev}.geo_map(:,:,iold), msh_new.geo_map);
        msh_lev{lev}.geo_map_jac(:,:,:,indices) = cat (4, hmsh.msh_lev{lev}.geo_map_jac(:,:,:,iold), msh_new.geo_map_jac);
        msh_lev{lev}.geo_map_der2(:,:,:,:,indices) = cat (5, hmsh.msh_lev{lev}.geo_map_der2(:,:,:,:,iold), msh_new.geo_map_der2);
        msh_lev{lev}.jacdet(:,indices) = [hmsh.msh_lev{lev}.jacdet(:,iold), msh_new.jacdet];
        msh_lev{lev}.element_size(:,indices) = [hmsh.msh_lev{lev}.element_size(:,iold), msh_new.element_size];
    end
end
end

function  neighbors  = get_neighbors_boundaries( hmsh, lev, ind_x, ind_y, ind_z, n_rb )
% x direction
neighbors_x_dir = linspace(ind_x-n_rb, ind_x+n_rb, 2*n_rb+1);
neighbors_x = zeros(size(neighbors_x_dir));
for i=1:numel(neighbors_x_dir)
    neighbors_x(i) =  neighbors_x_dir(i);
end
neighbors_x = unique(neighbors_x(neighbors_x>0));
% y direction
neighbors_y = [];
if numel(hmsh.mesh_of_level(lev).nel_dir) > 1
    neighbours_y_dir = linspace(ind_y-n_rb, ind_y+n_rb, 2*n_rb+1);
    neighbors_y = zeros(size(neighbours_y_dir));
    for i=1:numel(neighbours_y_dir)
        if(neighbours_y_dir(i)<=hmsh.mesh_of_level(lev).nel_dir(2)  && neighbours_y_dir(i)>0)
            neighbors_y(i) = sub2ind(hmsh.mesh_of_level(lev).nel_dir, ind_x, neighbours_y_dir(i), ind_z);
        end
    end
end
neighbors_y = unique(neighbors_y(neighbors_y>0));
neighbors = unique([neighbors_x neighbors_y]);
end


