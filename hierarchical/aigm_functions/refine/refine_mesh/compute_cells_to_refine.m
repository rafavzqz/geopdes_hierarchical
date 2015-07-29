function [ME, ind] = compute_cells_to_refine(hspace, hmsh, MF)
%
% function [ME, ind] = compute_cells_to_refine(hmsh, MF)
%
% This function computes the indices of cells that have to be splitted when
% marking for refinement the functions in MF.
%
% Input:    hspace:
%           hmsh:
%           MF: (cell array), where MF{lev} is a matrix whose rows are the indices
%                   of marked functions of level lev, for lev = 1:hspace.nlevels
%
% Output:   ME: (cell array),
%                   where ME{lev} is a matrix whose rows are the indices
%                   of marked elements of level lev, for lev = 1:hmsh.nlevels
%           ind: (cell array),
%                   where indices{lev} satisfies ME{lev} = hmsh.active{lev}(indices{lev}).
%
% This function uses sp_get_cells
%

ME = cell(hmsh.nlevels,1);
ind = cell(hmsh.nlevels,1);

for lev = 1:hmsh.nlevels
    ME{lev} = zeros(0,1);
    ind{lev} = zeros(0,1);
end

for lev = 1:hspace.nlevels
    if ~isempty(MF{lev})
        ME{lev} = sp_get_cells(hspace.space_of_level(lev), hmsh.mesh_of_level(lev), MF{lev});
        % Remove the nonactive cells from ME{lev}
        [ME{lev}, dummy, ind{lev}] = intersect(ME{lev}, hmsh.active{lev}, 'rows'); % ME{lev} = hmsh.active{lev}(ind{lev},:)
        % Now, ME{lev} contains the cells that have to be splitted
    end
end