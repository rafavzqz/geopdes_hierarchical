function [marked, nmarked] = adaptivity_mark(est, hmsh, hspace, adaptivity_data)
%
% function [marked, nmarked] = adaptivity_mark(est, hmsh, hspace, adaptivity_data)
%
% This function marks cells or basis functions for refinement according to 
% the marking strategy in adaptivity_data.mark_strategy. Possible strategies are
%           GR: global (uniform) refinement,  
%           MS: maximum strategy,  
%           GERS: guaranteed error reduction strategy (D\"orfler's)
%
% OUTPUT: marked{lev}: indices of marked cells (or functions) of level lev
%         nmarked: total number of marked cells (or functions)
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ATENCION: Despues describire mejor esta funcion y limpiare un poco mas el
% codigo
%

max_est = max(est);

switch adaptivity_data.flag
    case 'elements', disp('marking elements for refinement')
    case 'functions', disp('marking basis functions for refinement')
    otherwise, disp('MARK: Error'), return,
end

aux_marked = zeros(size(est));

switch adaptivity_data.mark_strategy
case 'GR'
  aux_marked = ones(size(est));
case 'MS'
  aux_marked(est > adaptivity_data.mark_param * max_est) = 1;
case 'GERS'
  est_sum2 = sum(est.^2);
  est_sum2_marked = 0;
  threshold = (1 - adaptivity_data.mark_param)^2 * est_sum2;
  gamma = 1;
  while (est_sum2_marked < threshold)
    gamma = gamma - 0.1;
    ff = find(est > gamma * max_est);
    aux_marked(ff) = 1;
    est_sum2_marked = sum((est(ff)).^2);
  end
end

marked_list = find(aux_marked);
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