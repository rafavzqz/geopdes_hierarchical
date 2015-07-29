function [hmsh, hspace] = refine (hmsh, hspace, marked, flag, boundary)
%
% function [hmsh, hspace] = refine (hmsh, hspace, marked, flag, boundary)
%
% This function updates the structures hmsh and hspace
% when enlarging the underlying subdomains with some marked functions or elements.
%
% Input:        hmsh: struct for current the hierarchical mesh
%               hspace: struct for current the hierarchical space
%               marked:  (cell array),
%                   where marked{lev} is a matrix whose rows are the indices
%                   of marked functions or elements of level lev, for lev =
%                   1,2,...
%               flag: 'functions' or 'elements'
%               boundary: true or false, default: true. (Fill the
%                   information for the boundaries of the mesh and space).
%
% Output:   hmsh: struct for the new hierarchical mesh after refinement
%           hspace: struct for the new hierarchical space after refinement
%
%

if nargin == 4
    boundary = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFINE MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

refine_mesh_time = tic;
disp('Refining mesh:')
switch flag
    case 'functions',
    [marked_elements, indices] = compute_cells_to_refine(hspace, hmsh, marked);
    case 'elements',
    marked_elements = marked;
    indices = [];
end
[hmsh, new_cells] = refine_hierarchical_mesh(hmsh, marked_elements, indices, boundary);
tempo = toc(refine_mesh_time);
fprintf('refine: Number of current active cells: %d (%f seconds)\n', hmsh.nel, tempo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFINE SPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

refine_space_time = tic;
disp('Updating space:')
hspace = refine_hierarchical_space(hmsh, hspace, marked, flag, new_cells, boundary);
tempo = toc(refine_space_time);
fprintf('Number of current dofs: %d (%f seconds)\n', hspace.ndof, tempo);
