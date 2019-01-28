% ELEMENTS_TO_BOXES: retrieve G+smo box data (see G+smo documentation) from
% elements indices (local to their corresponding levels).
%
% function boxes = elements_to_boxes (elements, hmsh)
%
% INPUT
%   elements: cell of size nlevels x 1 containing the elements of interest
%   hmsh: hierarchical_mesh object (see hierarchical_mesh.m)
% 
% OUTPUT
%   boxes: array of size 2*ndim+1 representing the given elements (see 
%       G+smo documentation of refine method)
%
% Copyright (C) 2018 Ondine Chanon

function boxes = elements_to_boxes (elements, hmsh)

elements = reshape (elements, [], 1);
cs = [0; cumsum(cellfun ('length', elements))];
boxes = zeros (cs(end), 2 * hmsh.ndim + 1);

for lev = 1:hmsh.nlevels
    if ~isempty (elements{lev})
        I = cell (hmsh.ndim,1);
        nel_dir_lev = hmsh.mesh_of_level(lev).nel_dir;
        [I{:}] = ind2sub (nel_dir_lev, elements{lev});
        I = cellfun (@(C) C-1, I, 'UniformOutput', false);
        I = reshape (cell2mat(I)', [], hmsh.ndim);
        I = [I, I+1];
        boxes(cs(lev)+1:cs(lev+1),:) = [lev * ones(size(I,1), 1), I];
    end
end

boxes = reshape (boxes', 1, []);

end

