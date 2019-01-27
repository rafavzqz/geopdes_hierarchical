% BOXES_TO_ELEMENTS: retrieve the elements indices (local to their
% corresponding levels) from G+smo box data.
%
% function elements = boxes_to_elements (lowpts, uppts, lev_boxes, nsub, nel_dirL1)
%
% INPUT
%   lowpts: [n x ndim] left bottom corners of boxes, with n = number of
%       boxes. The corners are described by their parametric index 
%       coordinates starting from 0.
%	uppts:	[n x ndim] right upper corners of boxes.
%	lev_boxes: [n x 1] corresponding levels in which the boxes lie
%       (starting from 1).
%   nsub: [ndim x 1] number of sub-elements to be added at each refinement 
%       step in each parametric direction (2 for dyadic).
%   nel_dirL1: [ndim x 1] number of elements in each parametric direction 
%       at level 1.
%
%   Note that the three first arguments, lowpts, uppts, lev_boxes correspond
%   to the outputs of the method getBoxes of a gsTHBSplineBasis (see G+smo 
%   documentation). 
% 
% OUTPUT
%   elements: cell of size [nlevels x 1] containing the local indices of 
%       the elements corresponding to the given boxes at their 
%       corresponding level
%
% Copyright (C) 2018 Ondine Chanon

function elements = boxes_to_elements (lowpts, uppts, lev_boxes, nsub, nel_dirL1)

if (~isequal (size (lowpts), size (uppts)) || length (lev_boxes) ~= size (lowpts,1) || ...
        length (nsub) ~= length (nel_dirL1) || length (nsub) ~= size (lowpts,2))
    error('Wrong input size')
end

maxLev = max (lev_boxes);
elements = cell (maxLev, 1);
nsub = reshape (nsub, [], 1);
nel_dirL1 = reshape (nel_dirL1, 1, []);
nel_dir_allL = nel_dirL1 .* (repmat(nsub, 1, maxLev) .^ ((0:maxLev-1)))';
nel_dir_allL = [nel_dir_allL ones(size(nel_dir_allL,1),1)];

nbox = length(lev_boxes);
lowpts_lev = lowpts ./ (repmat(nsub, 1, nbox)' .^ (maxLev - lev_boxes)) + 1; 
uppts_lev = uppts ./ (repmat(nsub, 1, nbox)' .^ (maxLev - lev_boxes));

for box = 1:nbox % loop over all the boxes
    lev = lev_boxes(box);
    index = cell (size (nel_dirL1));
    for ii = 1:length (nel_dirL1)
        index{ii} = lowpts_lev(box, ii):uppts_lev(box, ii);
    end
    all_idx = cell (size (nel_dirL1));
    [all_idx{:}] = ndgrid (index{:});
    all_idx = cellfun (@(C) reshape(C, 1, [])', all_idx, 'UniformOutput', false);

    elements{lev} = [elements{lev}; sub2ind(nel_dir_allL(lev,:), all_idx{:})]; 
end

elements = cellfun(@(C)sort(C), elements, 'UniformOutput', false);

end
