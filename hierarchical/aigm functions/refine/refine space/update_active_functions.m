function [A, W, R] = update_active_functions(A, W, R, E, RE, new_cells, I, degree, nel_dir, proj, ndof_dir, space_type)
%
% function [A, W, R] = update_active_functions(A, W, R, E, RE, new_cells, I, degree, nel_dir, proj, ndof_dir, space_type)
%
% Update active functions
%
% Input:    lev: (scalar) current level
%           A: (matrix) global indices of active functions of level lev (one row per function)
%           W: (column) coefficients for the partition-of-unity (one row per function)
%           R: (matrix) global indices of removed functions of level lev (one row per function)
%           EE: (matrix) global indices of active cells of level lev + 1 (one row per cell)
%           RREE: (matrix) global indices of refined cells of level lev + 1 (one row per cell)
%           new_cells: (matrix) global indices of the new cells of level lev + 1 after mesh refinement (one row per cell)
%           AA: (matrix) global indices of active functions of level lev +1 (one row per function)
%           WW: (column) coefficients for the partition-of-unity (one row per function)
%           RR: (matrix) global indices of removed functions of level lev + 1 (one row per function)
%           I: global indices of active functions to be refined
%           degree: array containing the polynomial degree in each coordinate direction
%           nelem_lev: array containing the total number of elements in the cartesian
% grid in this level in each coordinate direction
%
% Output:   [A, W, R, AA, WW, RR] updated
%
%
% This function uses:       split_basis
%                           get_cells
%

% Atencion: Despues podriamos trabajar directamente con un indice numerico en
% vez de indice tensorial para almacenar y modicar la info de A, R, y M

if size(A{1},2)~=size(I{1},2)
    disp('Error: Bad call to update_active_functions');
    return,
end

dim = size(A{1},2);

nlevels = numel(A);

if ~isempty(I{nlevels}) % if a new level is going to be activated
    A{nlevels+1,1} = zeros(0,dim);
    R{nlevels+1,1} = zeros(0,dim);
    W{nlevels+1,1} = zeros(0,1);
end

for lev = 1:nlevels
    I{lev}= union(I{lev}, intersect(R{lev}, A{lev}, 'rows'),'rows');
    if ~isempty(I{lev})
        [uno,indA] = ismember(I{lev}, A{lev},'rows');
        if any(uno~=1)
            disp('ERROR: update_active_functions: Some nonactive functions were selected');
        end
        
        coefficients = split_basis(proj(lev,:));
        nfunctions = size(I{lev},1);
        % Remove the corresponding functions from the active functions of level lev
        A{lev}(indA,:) = [];
        w = W{lev}(indA);
        W{lev}(indA) = [];
        R{lev} = union(R{lev}, I{lev}, 'rows'); % no funciona cuando R es vacio
        %R{lev} = vertcat(R{lev}, I{lev});
        %R{lev} = sortrows(R{lev});
        for i = 1: nfunctions
            %[Ichildren,c] = split_fun(I(i,:), lev);
            % Mejorar lo siguiente
            switch dim
                case 1,ind = I{lev}(i,1);
                case 2, ind = sub2ind(ndof_dir{lev}, I{lev}(i,1), I{lev}(i,2));
                case 3,ind = sub2ind(ndof_dir{lev}, I{lev}(i,1), I{lev}(i,2), I{lev}(i,3));
            end
            c = coefficients(:,ind);
            II = find(c);
            %c = kron(proj{2}(:,I(i,2)), proj{1}(:,I(i,1)));
            %II = find(abs(c) > 1e-10);
            
            c = c(II);
            c = full(c);
            % Mejorar lo siguiente
            switch dim
                case 1,Ichildren = ind2sub(ndof_dir{lev+1},II);
                    Ichildren = Ichildren(:);
                case 2, [a1,a2] = ind2sub(ndof_dir{lev+1},II);
                    Ichildren = [a1(:) a2(:)];
                case 3, [a1,a2,a3] = ind2sub(ndof_dir{lev+1},II);
                    Ichildren = [a1(:) a2(:) a3(:)];
            end
            % Update R{lev+1}: Computation of functions to be added to R{lev+1}
            if ~isempty(A{lev+1})
                Ichildren_nonactive = setdiff(Ichildren,A{lev+1},'rows');
            else
                Ichildren_nonactive = Ichildren;
            end
            if ~isempty(Ichildren_nonactive)
                if ~isempty(R{lev+1})
                    II = setdiff(Ichildren_nonactive,R{lev+1},'rows');
                else
                    II = Ichildren_nonactive;
                end
                if ~isempty(II)
                    nfun = size(II,1);
                    flag = zeros(1,nfun);
                    for ii = 1: nfun
                        flag(ii) = isempty(intersect(get_cells(II(ii,:), degree, nel_dir{lev+1}), E{lev+1},'rows'));
                    end
                    II = II(flag==1,:);
                    R{lev+1} = vertcat(R{lev+1}, II);
                end
                nchildren_nonactive = size(Ichildren_nonactive,1);
                A{lev+1} = vertcat(A{lev+1}, Ichildren_nonactive);
                W{lev+1} = vertcat(W{lev+1}, zeros(nchildren_nonactive,1));
            end
            [unos, indices] = ismember(Ichildren, A{lev+1}, 'rows');
            W{lev+1}(indices) = W{lev+1}(indices) + w(i)*c;
        end % for i = 1: nfunctions
        
        if space_type
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Now, we activate functions of level lev+1
            %% that are not children of any removed function of level lev
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            new_possible_active_fun = zeros(0,dim);
            number_of_new_cells = size(new_cells{lev+1},1);
            for el = 1:number_of_new_cells
                new_possible_active_fun = union(new_possible_active_fun,get_basis_functions(new_cells{lev+1}(el,:),degree, nel_dir{lev+1}),'rows');
            end
            new_possible_active_fun = setdiff(new_possible_active_fun,A{lev+1},'rows');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % We activate the new functions of level lev+1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            number_possible_new_functions = size(new_possible_active_fun,1);
            
            for fun = 1:number_possible_new_functions
                elem = get_cells(new_possible_active_fun(fun,:),degree, nel_dir{lev+1});
                if all( ismember(elem,E{lev+1},'rows') | ismember(elem,RE{lev+1},'rows') )
                    A{lev+1} = vertcat(A{lev+1}, new_possible_active_fun(fun,:));
                    W{lev+1} = vertcat(W{lev+1}, 0);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    if isempty(W{lev})
        A{lev} = zeros(0,dim);
        W{lev} = zeros(0,1);
    end
end