function [new_marked] = hrefine_rec(hmsh, hspace, Q_ind, lev_Q, m)
%REFINE_RECURSIVE algorithm
%hmsh: hierarchical mesh to be refined
%hspace: corresponding hierarchical space to be refined
%Q_ind: indices of the cell(s) to be refined
%lev_Q: level of the cell(s) Q_ind
%m: class of admissibility

new_marked = cell (1,hmsh.nlevels);
lev_i=lev_Q-m+1;
neighbours = neigh (hmsh,hspace, Q_ind, lev_Q, m);
if numel(neighbours)>0
    new_new_marked = hrefine_rec (hmsh, hspace, neighbours, lev_i, m);
    for k=1:lev_i
        new_marked{k} = union (new_marked{k}, new_new_marked{k});
    end
end

if (numel(intersect(Q_ind,hmsh.active{lev_Q})) > 0)
    new_marked{lev_Q} = union (new_marked{lev_Q}, Q_ind);
end

end