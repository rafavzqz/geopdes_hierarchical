function [hmsh, hspace] = refine (hmsh, hspace, marked, flag, flag_full_basis, boundary)
%
% function [hmsh, hspace] = refine (hmsh, hspace, marked, flag, flag_whole_basis, boundary)
%
% This function updates the structures hmsh and hspace
% when enlarging the underlying subdomains with some marked functions or elements.
%
% Input:        hmsh: struct for current the hierarchical mesh
%               hspace: struct for current the hierarchical space
%               marked:  (1 x msh.nlevels cell array),
%                   where marked{lev} is a matrix whose rows are the tensor-product indices
%                   of marked functions or elements of level lev, for lev =
%                   1:hmsh.nlevels
%               flag: 'functions' or 'elements'
%               flag_full_basis: true or false
%               boundary: true or false, default: true. (Fill the
%                   information for the boundaries of the mesh and space).
%
% Output:   hmsh: struct for the new hierarchical mesh after refinement
%           hspace: struct for the new hierarchical space after refinement
%
%

if nargin == 5
    boundary = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REFINE MESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

refine_mesh_time = tic;
disp('Refining mesh:')
switch flag
    case 'functions',
    [marked_elements, indices] = compute_cells_to_refine(hmsh, marked, hspace.degree);
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
switch flag
    case 'functions',
        functions_to_remove = marked;
    case 'elements',
        functions_to_remove = compute_functions_to_remove(hmsh, hspace, marked, flag);
end
hspace = refine_hierarchical_space(hmsh, hspace, functions_to_remove, new_cells, flag_full_basis, boundary);
tempo = toc(refine_space_time);
fprintf('Number of current dofs: %d (%f seconds)\n', hspace.ndof, tempo);
