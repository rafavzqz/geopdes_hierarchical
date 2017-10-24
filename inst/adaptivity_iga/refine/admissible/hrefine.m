function [marked] = hrefine( hmsh, hspace, toberef, m)
%REFINE algorithm
%hmsh: hierarchical mesh to be refined
%hspace: corresponding hierarchical space to be refined
%toberef: cell array containing the cells to be refined (one vector for each level)
%m: class of admissibility
if m<2
    error('class of admissibility m<2')
end

marked=cell(1,hmsh.nlevels);
for l=1:hmsh.nlevels
    Q_ind=intersect(toberef{l},hmsh.active{l});
    if numel(Q_ind)>0
        new_marked=hrefine_rec(hmsh, hspace, Q_ind, l, m);
        for k=1:l
            marked{k}=union(marked{k}, new_marked{k});
        end
    end
end
end