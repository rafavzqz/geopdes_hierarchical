function I = compute_functions_to_remove(hmsh, hspace, marked, flag)
%
% function I = compute_functions_to_remove(hmsh, hspace, marked, flag)
%
% Computation of active functions to be removed in each level.
% ATENCION: Voy a limpiar y editar esta funcion
%
% Input:
%           hmsh:
%           hspace:
%           marked{lev}: (matrix) global indices of marked functions or
%           elements of level lev, for lev = 1:nlevels
%           flag: 'functions' or 'elements', default: 'functions'. flag
%           indicates if marked is a set of functions or elements.
%
% Output:   I{lev}: global indices of active functions of level lev to be
% removed, for lev = 1:nlevels
%
%
% This function uses:   get_neighbors
%                       get_cells
%                       get_basis_functions
%

% We fill E with the information of active cells
Ne = cumsum([0; hmsh.nel_per_level(:)]);
E = cell(hmsh.nlevels+1,1);
% El siguiente loop se puede evitar usando mat2cell
for lev = 1:hmsh.nlevels
        ind_e = (Ne(lev)+1):Ne(lev+1);
        E{lev} = hmsh.globnum_active(ind_e, 2:end);
end

% We fill A with the information of active functions
Nf = cumsum([0; hspace.ndof_per_level(:)]);
A = cell(hspace.nlevels+1,1);
% El siguiente loop se puede evitar usando mat2cell
for lev = 1:hspace.nlevels
        ind_f = (Nf(lev)+1):Nf(lev+1);
        A{lev} = hspace.globnum_active(ind_f, 2:end);
end


I = cell(1, hspace.nlevels);
for lev = 1:hspace.nlevels
    M = marked{lev};
    nmarked = size(M,1);
    switch flag,
        case 'functions',
            % Computation of functions which possibly have to be removed
            %%%%%%%%%%%%%%%%%%%%%%
            % Mejorar el siguiente bloque
            %%%%%%%%%%%%%%%%%%%%%%
            I{lev} = [];
            for i = 1: nmarked
                aux = get_neighbors(M(i,:), hspace.degree, hmsh.mesh_of_level(lev).nel_dir);
                I{lev} = [I{lev}; aux];
            end
            I{lev} = unique(I{lev},'rows');
            %%%%%%%%%%%%%%%%%%%%%%
            % Remove from I the nonactive functions and the active functions already
            % selected for refinement
            if ~isempty(M)
                I{lev} = setdiff(intersect(I{lev}, A{lev}, 'rows'), M, 'rows');
            else
                I{lev} = intersect(I{lev}, A{lev}, 'rows');
            end
        case 'elements',
            % Computation of functions which possibly have to be removed
            %%%%%%%%%%%%%%%%%%%%%%
            % Mejorar el siguiente bloque
            %%%%%%%%%%%%%%%%%%%%%%
            I{lev} = [];
            for i = 1:nmarked
                aux = get_basis_functions(M(i,:), hspace.degree, hmsh.mesh_of_level(lev).nel_dir);
                I{lev} = [I{lev}; aux];
            end
            I{lev} = unique(I{lev},'rows');
            %%%%%%%%%%%%%%%%%%%%%%
            % Remove from I the nonactive functions
            I{lev} = intersect(I{lev}, A{lev}, 'rows');
    end
    
    % Computation of functions which in fact have to be removed
    nfunctions = size(I{lev},1);
    flag_ell = zeros(1,nfunctions);
    for i = 1: nfunctions
        flag_ell(i) = isempty(intersect(get_cells(I{lev}(i,:), hspace.degree, hmsh.mesh_of_level(lev).nel_dir), E{lev}, 'rows'));
    end
    I{lev} = I{lev}(flag_ell == 1,:);
    
    if strcmp(flag,'functions')
        I{lev} = [M; I{lev}];
    end
    
end