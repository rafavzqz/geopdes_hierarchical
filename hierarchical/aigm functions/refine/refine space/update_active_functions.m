function hspace = update_active_functions(hspace, hmsh, new_cells, I)
%
% function hspace = update_active_functions(hspace, hmsh, new_cells, I)
%
% This function updates the active dofs (hspace.active and hspace.globnum_active), their coefficients (hspace.coeff) and deactivated dofs (hspace.deactivated) in each level when
% refining the functions in I. This function also updates hspace.nlevels, hspace.ndof and hspace.ndof_per_level
% ATENCION: Voy a mejorar un poco la performance de esta funcion
%
% Input:    hspace:
%           hmsh:
%           new_cells: see refine_hierarchical_mesh
%           I{lev}: indices of active functions of level lev to be deactivated
%
% Output:   hspace:
%
%
% This function uses:       split_basis
%                           sp_get_cells
%                           sp_get_basis_functions
%

% Mejorar este chequeo
if size(hspace.active{1},2)~=size(I{1},2)
    disp('Error: Bad call to update_active_functions');
    return,
end

Nf = cumsum([0 hspace.ndof_per_level]);
W = cell(hspace.nlevels+1,1);
active = cell(hspace.nlevels+1,1);
deactivated = cell(hspace.nlevels+1,1);
% El siguiente loop seguramente se puede evitar usando mat2cell
for lev = 1:hspace.nlevels
    ind_f = (Nf(lev)+1):Nf(lev+1);
    active{lev} = hspace.active{lev};
    deactivated{lev} = hspace.deactivated{lev};
    W{lev} = hspace.coeff(ind_f);
    W{lev} = W{lev}(:);
end

active{hspace.nlevels+1,1} = zeros(0,1);
deactivated{hspace.nlevels+1,1} = zeros(0,1);
W{hspace.nlevels+1,1} = zeros(0,1);

if hspace.nlevels == hmsh.nlevels
    max_lev = hspace.nlevels - 1;
else
    max_lev = hspace.nlevels;
end

for lev = 1:max_lev
    I{lev}= union(I{lev}, intersect(deactivated{lev}, active{lev}, 'rows'),'rows');
    if ~isempty(I{lev})
        [uno,indA] = ismember(I{lev}, active{lev},'rows');
        if any(uno~=1)
            disp('ERROR: update_active_functions: Some nonactive functions were selected');
        end
        
        coefficients = split_basis(hspace.Proj(lev,:));
        nfunctions = size(I{lev},1);
        % Remove the corresponding functions from the active functions of level lev
        active{lev}(indA) = [];
        w = W{lev}(indA);
        W{lev}(indA) = [];
        deactivated{lev} = union(deactivated{lev}, I{lev}, 'rows');
        % Vectorizar lo que se pueda en el siguiente loop, en particular
        % para hacer un solo llamado a sp_get_cells
        for i = 1: nfunctions
            ind = I{lev}(i,:);
            c = coefficients(:,ind);
            II = find(c);
            c = c(II);
            c = full(c);
            Ichildren = II;
            % Update deactivated{lev+1}: Computation of functions to be
            % added to deactivated{lev+1}
            Ichildren_nonactive = setdiff(Ichildren,active{lev+1},'rows');
            if ~isempty(Ichildren_nonactive)
                II = setdiff(Ichildren_nonactive,deactivated{lev+1},'rows');
                if ~isempty(II)
                    nfun = size(II,1);
                    flag = zeros(1,nfun);
                    [dummy, cells_per_fun] = sp_get_cells(hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), II);
                    for ii = 1: nfun
                        flag(ii) = isempty(intersect(cells_per_fun{ii}, hmsh.active{lev+1},'rows'));
                    end
                    II = II(flag==1,:);
                    deactivated{lev+1} = vertcat(deactivated{lev+1}, II);
                end
                nchildren_nonactive = size(Ichildren_nonactive,1);
                active{lev+1} = vertcat(active{lev+1}, Ichildren_nonactive);
                W{lev+1} = vertcat(W{lev+1}, zeros(nchildren_nonactive,1));
            end
            [unos, indices] = ismember(Ichildren, active{lev+1}, 'rows');
            W{lev+1}(indices) = W{lev+1}(indices) + w(i)*c;
        end % for i = 1: nfunctions
    end
    if (hspace.type && ~isempty(new_cells{lev+1}))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Now, we activate functions of level lev+1
        %% that are not children of any removed function of level lev
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        new_possible_active_fun = sp_get_basis_functions(hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1),new_cells{lev+1});
        new_possible_active_fun = setdiff(new_possible_active_fun,active{lev+1},'rows');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % We activate the new functions of level lev+1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        number_possible_new_functions = size(new_possible_active_fun,1);
        
        [dummy, elem] = sp_get_cells(hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1),new_possible_active_fun);
        
        for fun = 1:number_possible_new_functions
            if all( ismember(elem{fun},hmsh.active{lev+1},'rows') | ismember(elem{fun},hmsh.deactivated{lev+1},'rows') )
                active{lev+1} = vertcat(active{lev+1}, new_possible_active_fun(fun,:));
                W{lev+1} = vertcat(W{lev+1}, 0);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if isempty(W{lev})
        active{lev} = zeros(0,1);
        W{lev} = zeros(0,1);
    end
end

if ~isempty(active{hspace.nlevels + 1})
    hspace.nlevels = hspace.nlevels + 1;
end

hspace.active = cell(hspace.nlevels,1);
hspace.deactivated = cell(hspace.nlevels,1);
for lev = 1: hspace.nlevels
    hspace.active{lev} = active{lev};
    hspace.deactivated{lev} = deactivated{lev};
end

hspace.coeff = cell2mat(W);


% We update hspace.ndof_per_level
ndof_per_level = zeros(1, hspace.nlevels);
%aux_ind_fun = [];
for lev = 1:hspace.nlevels
    ndof_per_level(lev) = size(hspace.active{lev},1);
    %aux_ind_fun = vertcat(aux_ind_fun, lev*ones(ndof_per_level(lev),1));
end

hspace.ndof_per_level = ndof_per_level;

% We update hspace.ndof
hspace.ndof = sum(hspace.ndof_per_level);

% We update hspace.globnum_active
hspace.globnum_active = zeros(hspace.ndof,hmsh.ndim+1);
Nf = cumsum([0; hspace.ndof_per_level(:)]);
for lev = 1:hspace.nlevels
    ind_f = (Nf(lev)+1):Nf(lev+1);
    if ~isempty(ind_f)
        hspace.globnum_active(ind_f, 1) = lev;
        globnum = cell(1,hmsh.ndim);
        [globnum{:}] = ind2sub(hspace.space_of_level(lev).ndof_dir,hspace.active{lev}(:));
        hspace.globnum_active(ind_f, 2:end) = cell2mat(globnum);
    end
end
