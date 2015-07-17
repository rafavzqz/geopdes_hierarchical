function [A, W, R, AA, WW, RR] = update_active_functions(lev, A, W, R, EE, RREE, new_cells, AA, WW, RR, I, degree, nelem_lev, proj, current_ndof_dir, new_ndof_dir,flag_whole_basis)
%
% function [A, W, R, AA, WW, RR] = update_active_functions(lev, A, W, R, EE, RREE, new_cells, AA, WW, RR, I, degree, nelem_lev, proj)
%
% Update active functions (remove functions of level lev and add active functions of level lev + 1).
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

I= [I; intersect(R, A, 'rows')];
I = unique(I, 'rows');
[uno,indA] = ismember(I, A,'rows');
if any(uno~=1)
    disp('ERROR: get_functions_to_remove: Some nonactive functions were selected');
end

dim = numel(degree);
coefficients = split_basis(proj);
nfunctions = size(I,1);
 % Remove the corresponding functions from the active functions of level lev
A(indA,:) = [];
w = W(indA);
W(indA) = [];
% R = union(R, I, 'rows'); % no funciona cuando R es vacio
R = vertcat(R, I);
R = sortrows(R); 
% n1 = size(AA,1); % initial number of active (and in AA\cap RR) functions of level lev + 1
for i = 1: nfunctions
    %[Ichildren,c] = split_fun(I(i,:), lev);
    % Mejorar lo siguiente
    switch dim
        case 1,ind = I(i,1);
        case 2, ind = sub2ind(current_ndof_dir, I(i,1), I(i,2));
        case 3,ind = sub2ind(current_ndof_dir, I(i,1), I(i,2), I(i,3));
    end
    c = coefficients(:,ind);
    II = find(c);
    %c = kron(proj{2}(:,I(i,2)), proj{1}(:,I(i,1)));
    %II = find(abs(c) > 1e-10);
    
    c = c(II);
    c = full(c);
    % Mejorar lo siguiente
    switch dim
        case 1,Ichildren = ind2sub(new_ndof_dir,II);
            Ichildren = Ichildren(:);
        case 2, [a1,a2] = ind2sub(new_ndof_dir,II);
            Ichildren = [a1(:) a2(:)];
        case 3, [a1,a2,a3] = ind2sub(new_ndof_dir,II);
            Ichildren = [a1(:) a2(:) a3(:)];
    end
    % Update RR: Computation of functions to be added to RR
    if ~isempty(AA)
        Ichildren_nonactive = setdiff(Ichildren,AA,'rows');
    else
        Ichildren_nonactive = Ichildren;
    end
    if ~isempty(Ichildren_nonactive)
        if ~isempty(RR)
            II = setdiff(Ichildren_nonactive,RR,'rows');
        else
            II = Ichildren_nonactive;
        end
        if ~isempty(II)
            nfun = size(II,1);
            flag = zeros(1,nfun);
            for ii = 1: nfun
                flag(ii) = isempty(intersect(get_cells(II(ii,:), degree, 2*nelem_lev), EE,'rows'));
            end
            II = II(flag==1,:);
            RR = vertcat(RR, II);
        end
        nchildren_nonactive = size(Ichildren_nonactive,1);
        AA = vertcat(AA, Ichildren_nonactive);
        WW = vertcat(WW, zeros(nchildren_nonactive,1));
    end
    [unos, indices] = ismember(Ichildren, AA, 'rows');
    WW(indices) = WW(indices) + w(i)*c;
end % for i = 1: nfunctions
% n2 =  size(AA,1); % final number of active (and in AA\cap RR) functions of level lev + 1

if flag_whole_basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now, we activate functions of level lev+1 
%% that are not children of any removed function of level lev 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

new_possible_active_fun = zeros(0,dim);
number_of_new_cells = size(new_cells,1);
for el = 1:number_of_new_cells
    new_possible_active_fun = union(new_possible_active_fun,get_basis_functions(new_cells(el,:),degree,2*nelem_lev),'rows');
end
new_possible_active_fun = setdiff(new_possible_active_fun,AA,'rows');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We activate the new functions of level lev+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
number_possible_new_functions = size(new_possible_active_fun,1);

for fun = 1:number_possible_new_functions 
    elem = get_cells(new_possible_active_fun(fun,:),degree, 2*nelem_lev);
    if all( ismember(elem,EE,'rows') | ismember(elem,RREE,'rows') )
        AA = vertcat(AA, new_possible_active_fun(fun,:));
        WW = vertcat(WW, 0);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

