function [marked, nmarked] = marking_for_coarsening (est, hmsh, hspace, adaptivity_data)
%
% function [marked, nmarked] = marking_for_coarsening (est, hmsh, hspace, adaptivity_data)
%

switch adaptivity_data.flag
    case 'elements' %, disp ('marking elements for refinement')
    case 'functions' %, disp ('marking basis functions for refinement')
    otherwise, error ('adaptivity_mark: Unknown option %s', adaptivity_data.flag)
end

aux_marked = zeros (size (est));

switch adaptivity_data.mark_strategy
    case 'GR'
        aux_marked = ones (size (est));
    case 'MS'
        n = round(adaptivity_data.mark_param_coarsening*numel(est));
        [~, ind] = sort(est);
        aux_marked(ind(1:n)) = 1;
end

marked_list = find (aux_marked);
nmarked = numel (marked_list);

marked = cell (hmsh.nlevels, 1);

switch (lower (adaptivity_data.flag))
    case 'elements'
        aux = cumsum ([0, hmsh.nel_per_level]);
        for lev = 1:hmsh.nlevels
            elems = aux(lev)+1:aux(lev+1);
            [~,ind,~] = intersect (elems, marked_list);
            marked{lev} = hmsh.active{lev}(ind);
        end
    case 'functions'
        aux = cumsum ([0, hspace.ndof_per_level]);
        for lev = 1:hmsh.nlevels
            funs = aux(lev)+1:aux(lev+1);
            [~,ind,~] = intersect (funs, marked_list);
            marked{lev} = hspace.active{lev}(ind);
        end
end

end