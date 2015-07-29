function I = compute_functions_to_deactivate(hmsh, hspace, M, flag)
%
% function I = compute_functions_to_deactivate(hmsh, hspace, M, flag)
%
% Computation of active functions to be removed in each level.
%
% Input:
%           hmsh:
%           hspace:
%           M{lev}: global indices of marked functions or
%           elements of level lev, for lev = 1,2,...
%           flag: 'functions' or 'elements'. flag
%           indicates if marked is a set of functions or elements.
%
% Output:   I{lev}: global indices of active functions of level lev to be
% removed, for lev = 1:hspace.nlevels
%

% This function uses:   sp_get_neighbors
%                       sp_get_cells
%                       sp_get_basis_functions
%

I = cell(hspace.nlevels,1);

for i = 1:hspace.nlevels
    I{i} = zeros(0,1);
end

for lev = 1:hspace.nlevels
    % Computation of functions which possibly have to be deactivated
    if ~isempty(M{lev})
        switch flag,
            case 'functions',
                I{lev} = sp_get_neighbors(hspace.space_of_level(lev), hmsh.mesh_of_level(lev), M{lev});
                % Remove from I{lev} the nonactive functions and the active functions already
                % selected for refinement
                I{lev} = setdiff(intersect(I{lev}, hspace.active{lev}, 'rows'), M{lev}, 'rows');
            case 'elements',
                I{lev} =sp_get_basis_functions(hspace.space_of_level(lev), hmsh.mesh_of_level(lev), M{lev});
                % Remove from I{lev} the nonactive functions
                I{lev} = intersect(I{lev}, hspace.active{lev}, 'rows');
        end
        % Computation of functions that in fact have to be deactivated
        [dummy, cells_per_fun] = sp_get_cells(hspace.space_of_level(lev), hmsh.mesh_of_level(lev), I{lev});
        nfunctions = size(I{lev},1);
        flag_ell = zeros(1,nfunctions);
        for i = 1: nfunctions
            flag_ell(i) = isempty(intersect(cells_per_fun{i}, hmsh.active{lev}, 'rows'));
        end
        I{lev} = I{lev}(flag_ell == 1,:);
        
        if strcmp(flag,'functions')
            I{lev} = union(I{lev}, M{lev},'rows');
        end
    end
end