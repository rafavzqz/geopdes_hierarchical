function [deact_marked, num] = active2deactivated_marking(marked, hmsh, hspace, adaptivity_data)
%
% function [deact_marked, num] = active2deactivated_marking(marked, hmsh, hspace, adaptivity_data)
%
% INPUT:
%
%   marked:  cell-array with the indices of marked cells (or functions) for each level, in the tensor-product setting
%   hmsh:   object representing the coarse hierarchical mesh (see hierarchical_mesh)
%   hspace: object representing the coarse space of hierarchical splines (see hierarchical_space)
%   adaptivity_data: a structure with the data for the adaptivity method.
%                    In particular, it contains the field:
%                     -'flag': elements or functions, according to marked
%
% OUTPUT:
%
%    deact_marked:  cell-array with the indices of marked cells (or functions) for each level to be reactivated, in the tensor-product setting
%    num         :  number of cells (or elements) to be reactivated
%
% This function computes, acording to adaptivity_data.flag, the deactivated entities to be possibly reactivated. 
% If flag = elements, all their children are (active and) marked. 
% If flag = functions, they have at least one child (active and) marked.  
%


deact_marked = cell (hmsh.nlevels, 1);
for lev = 1:hmsh.nlevels-1
    if ~isempty(marked{lev+1})
        switch (lower (adaptivity_data.flag))
            case 'elements'
                [parents, flag] = hmsh_get_parent (hmsh, lev+1, marked{lev+1});
                if flag ~= 1
                    disp('Some nonactive elements were marked.')
                    return;
                end
                parents = intersect(parents, hmsh.deactivated{lev});
                for i = 1:numel(parents)
                    if all(ismember(hmsh_get_children (hmsh, lev, parents(i)), marked{lev+1}))
                        deact_marked{lev} = union(deact_marked{lev}, parents(i));
                    end
                end
            case 'functions'
                [parents, flag] = hspace_get_parents (hspace, lev+1, marked{lev+1});
                if flag ~= 1
                    disp('Some nonactive functions were marked.')
                    return;
                end
                parents = intersect(parents, hspace.deactivated{lev});
                for i = 1:numel(parents)
                    if any(ismember(hspace_get_children (hspace, lev, parents(i)), marked{lev+1}))
                        deact_marked{lev} = union(deact_marked{lev}, parents(i));
                    end
                end
        end
    end
    
end

num = numel(cell2mat(deact_marked'));