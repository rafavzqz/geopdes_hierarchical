function [M, ind] = compute_cells_to_refine(hmsh, M, degree)
%
% function [M, ind] = compute_cells_to_refine(hmsh, M, degree)
%
% This function computes the indices of cells that have to be splitted when
% marking for refinement the functions in M.
% ATENCION: Voy a modificar y limpiar este codigo cuando tenga la version
% general de get_cells
%
% Input:    hmsh:
%           M: (1 x msh.nlevels cell array), where M{lev} is a matrix whose rows are the tensor-product indices
%                   of marked functions of level lev, for lev = 1:hmsh.nlevels
%           degree: array containing the polynomial degree in each coordinate direction
%
% Output:   M: (1 x msh.nlevels cell array),
%                   where M{lev} is a matrix whose rows are the tensor-product indices
%                   of marked elements of level lev, for lev = 1:hmsh.nlevels
%           indices: (1 x msh.nlevels cell array),
%                   where indices{lev} satisfies M{lev} = hmsh.active{lev}(indices{lev}).
%
% This function uses get_cells
%

% We fill E with the information of active cells
Ne = cumsum([0; hmsh.nel_per_level(:)]);
E = cell(hmsh.nlevels+1,1);
% El siguiente loop se puede evitar usando mat2cell
for lev = 1:hmsh.nlevels
    if (hmsh.msh_lev{lev}.nel ~= 0)
        ind_e = (Ne(lev)+1):Ne(lev+1);
        E{lev} = hmsh.globnum_active(ind_e, 2:end);
    else
        E{lev} = zeros(0,hmsh.ndim);
    end
end

ind = cell(1,hmsh.nlevels);

for lev = 1:hmsh.nlevels
    if ~isempty(M{lev})
        nmarked = size(M{lev},1); % Number of marked functions of a fixed level
        nelem_lev = hmsh.mesh_of_level(lev).nel_dir;
        %%%%%%%%%%%%%%%%%%%%%%
        % Mejorar el siguiente bloque
        %%%%%%%%%%%%%%%%%%%%%%
        I = [];
        for i = 1: nmarked % For each marked function
            aux = get_cells(M{lev}(i,:), degree, nelem_lev);
            I = [I; aux];
        end
        I = unique(I,'rows');
        %%%%%%%%%%%%%%%%%%%%%%
        % Remove the nonactive cells from I
        [M{lev}, basura, ind{lev}] = intersect(I, E{lev}, 'rows'); % M{lev} = E{lev}(ind{lev},:)
        % Now, I contains the cells that have to be splitted
    end
end
