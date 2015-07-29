function [marked, nmarked] = mark(est, hmsh, hspace, adaptivity_data)
%
% function [marked, nmarked] = mark(est, hmsh, hspace, adaptivity_data)
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

nmarked = numel(marked_list);
switch adaptivity_data.flag
    case 'elements',
        globnum_marked = hmsh.globnum_active(marked_list,:);
        marked_sub = cell(hmsh.nlevels,1);
    case 'functions',
        globnum_marked = hspace.globnum_active(marked_list,:);
        marked_sub = cell(hspace.nlevels,1);
end

% Mejorar el siguiente procedimiento %%%%%%
for fila = 1:nmarked
    aux = globnum_marked(fila,1);
    marked_sub{aux} = [marked_sub{aux}; globnum_marked(fila,2:end)];
end

for lev = 1:numel(marked_sub)
    if isempty(marked_sub{lev})
        marked_sub{lev} = zeros(0,hmsh.ndim);
    end
end

%marked.sub = marked_sub;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

marked.ind = cell(size(marked_sub));

for lev = 1:numel(marked_sub)
    if isempty(marked_sub{lev})
        marked.ind{lev} = zeros(0,1);
    else
        switch adaptivity_data.flag
            case 'elements',
                siz = hmsh.mesh_of_level(lev).nel_dir;
            case 'functions',
                siz = hspace.space_of_level(lev).ndof_dir;
        end
        switch hmsh.ndim
            %case 1, marked_ind{lev} = marked_sub{lev}(:,1);
            case 2, marked.ind{lev} = sub2ind(siz, marked_sub{lev}(:,1), marked_sub{lev}(:,2));
            case 3, marked.ind{lev} = sub2ind(siz, marked_sub{lev}(:,1), marked_sub{lev}(:,2), marked_sub{lev}(:,3));
        end
    end
end

marked = marked.ind;