function [marked, nmarked] = mark(est, hmsh, hspace, mark_strategy, param,flag)
%
% function [marked, nmarked] = mark(est, hmsh, hspace, mark_strategy, param,flag)
%
% marks elements or basis functions for refinement according to 
% the marking strategy mark_strategy
% Possible strategies are
% GR: global (uniform) refinement,  
% MS: maximum strategy,  
% GERS: guaranteed error reduction strategy (D\"orfler's)

%nest = numel(est);
max_est = max(est);

switch flag
    case 'elements', disp('marking elements for refinement')
    case 'functions', disp('marking basis functions for refinement')
    otherwise, disp('MARK: Error'), return,
end

marked = zeros(size(est));

switch mark_strategy
case 'GR'
  marked = ones(size(est));
case 'MS'
  marked(est > param * max_est) = 1;
case 'GERS'
  est_sum2 = sum(est.^2);
  est_sum2_marked = 0;
  threshold = (1 - param)^2 * est_sum2;
  gamma = 1;
  while (est_sum2_marked < threshold)
    gamma = gamma - 0.1;
    ff = find(est > gamma * max_est);
    marked(ff) = 1;
    est_sum2_marked = sum((est(ff)).^2);
  end
end

marked_list = find(marked);

nmarked = numel(marked_list);
switch flag
    case 'elements',
        globnum_marked = hmsh.globnum_active(marked_list,:);
        marked = cell(hmsh.nlevels,1);
        ndim = hmsh.ndim;
    case 'functions',
        globnum_marked = hspace.globnum_active(marked_list,:);
        marked = cell(hspace.nlevels,1);
        ndim = hspace.ndim;
end

% Mejorar el siguiente procedimiento %%%%%%
for fila = 1:nmarked
    aux = globnum_marked(fila,1);
    marked{aux} = [marked{aux}; globnum_marked(fila,2:end)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:numel(marked)
    if isempty(marked{ii})
        marked{ii} = zeros(0,ndim);
    end
end




