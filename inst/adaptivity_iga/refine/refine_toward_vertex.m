function [ ref, nref_el ] = refine_toward_vertex( vertex, hmsh, hspace )
%REFINE_TOWARD_VERTEX Summary of this function goes here
%   Detailed explanation goes here

msh_finer_lev = hmsh.mesh_of_level(end);
sp_finer_lev = hspace.space_of_level(end);
aux_ref = zeros(sp_finer_lev.ndof,1);
local_vertex = msh_finer_lev.map(vertex);

% MARK ELEMENT OF THE FINER LEVEL CONATINING THE VERTEX
isInside = cell(hmsh.rdim,1);
% LOOP OVER DIRECTIONS
for idir = 1:hmsh.rdim
    isInside{idir} = zeros(1,msh_finer_lev.nel_dir(idir));
    % LOOP OVER ACTIVE DOFS OF THE FINER LEVEL
    for el = 1:msh_finer_lev.nel_dir(idir)
        degree_dir = sp_finer_lev.degree(idir);
        knots_dir = sp_finer_lev.knots{idir};
        % check if vertex is inside te knot-span in ith direction
        if (local_vertex(idir) > knots_dir(el+degree_dir) &&...
                local_vertex(idir) <= knots_dir(el+degree_dir+1))
            isInside{idir}(el) = 1;
        else
            isInside{idir}(el) = 0;
        end
    end % END ELEMENTS LOOP
end % END DIRECTION LOOP

index_ref = 1;

% LOOP OVER DIRECTIONS
for idir = 1:hmsh.rdim
    i = find( isInside{idir});
    if idir == 1
        index_ref = index_ref + (i-1);
    elseif idir == 2
        index_ref = index_ref + numel(isInside{1})*(i-1);
    else
        index_ref = index_ref + numel(isInside{1})*numel(isInside{2})*(i-1);
    end
        
end % END DIRECTION LOOP

% [cell_indices, ~] = sp_get_cells (sp_finer_lev, msh_finer_lev, index_ref);

aux_ref(index_ref) = 1;

ref_list = find (aux_ref);
nref_el = numel (ref_list);
ref = cell (hmsh.nlevels, 1);

aux = cumsum ([0, hmsh.nel_per_level]);

elems = aux(hmsh.nlevels)+1:aux(hmsh.nlevels+1); % active elements of the finer level
[~,ind,~] = intersect (elems, ref_list);
ref{end} = hmsh.active{end}(ind);

end

