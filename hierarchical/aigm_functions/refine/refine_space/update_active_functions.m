function hspace = update_active_functions (hspace, hmsh, new_cells, I)
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
% XXXX Change the help. Decide where to put this function (inside hspace_refine?)

% XXXX Mejorar este chequeo
if size(hspace.active{1},2)~=size(I{1},2)
    disp('ERROR: Bad call to update_active_functions');
    return,
end

Nf = cumsum ([0; hspace.ndof_per_level(:)]);
W = cell (hspace.nlevels+1,1);

active = hspace.active;
deactivated = hspace.deactivated;

% El siguiente loop seguramente se puede evitar usando mat2cell
for lev = 1:hspace.nlevels
  ind_f = (Nf(lev)+1):Nf(lev+1);
  W{lev} = hspace.coeff(ind_f);
  W{lev} = W{lev}(:);
end

active{hspace.nlevels+1,1} = zeros (0,1);
deactivated{hspace.nlevels+1,1} = zeros (0,1);
W{hspace.nlevels+1,1} = zeros (0,1);

if (hspace.nlevels == hmsh.nlevels)
    max_lev = hspace.nlevels - 1;
else
    max_lev = hspace.nlevels;
end

for lev = 1:max_lev
  I{lev} = union (I{lev}, intersect (deactivated{lev}, active{lev}));
  if (~isempty(I{lev}))
    [uno,indA] = ismember (I{lev}, active{lev});
    if (~all (uno))
      disp('ERROR: update_active_functions: Some nonactive functions were selected');
    end

    coefficients = 1;
    for idim = 1:hmsh.ndim
      coefficients = kron (hspace.Proj{lev,idim}, coefficients);
    end
    
    % Remove the corresponding functions from the active functions of level lev
    active{lev}(indA) = [];
    w = W{lev}(indA);
    W{lev}(indA) = [];
    deactivated{lev} = union (deactivated{lev}, I{lev});
    deactivated{lev} = deactivated{lev}(:);
    
    % Vectorizar lo que se pueda en el siguiente loop, en particular
    % para hacer un solo llamado a sp_get_cells
    for ifun = 1:numel (I{lev})
      ind = I{lev}(ifun,:);
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
      W{lev+1}(indices) = W{lev+1}(indices) + w(ifun)*c;
    end % for i = 1:numel(I{lev})
  end
  
% For the classical hierarchical space, we activate functions of level lev+1, 
%  that are not children of any removed function of level lev
% XXXX Change the way to use hspace.type
  if (hspace.type && ~isempty (new_cells{lev+1}))
        
    new_possible_active_fun = sp_get_basis_functions (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), new_cells{lev+1});
%     new_possible_active_fun = setdiff (new_possible_active_fun, active{lev+1}, 'rows');
% XXXXXX I THINK THIS CHANGE (the union) WAS NECESSARY. ASK EDUARDO IF IT IS OK
    new_possible_active_fun = setdiff (new_possible_active_fun, union (active{lev+1}, deactivated{lev+1}));

    [dummy, elem] = sp_get_cells (hspace.space_of_level(lev+1), hmsh.mesh_of_level(lev+1), new_possible_active_fun);

    new_functions = cellfun (@(x) all (ismember (x, union (hmsh.active{lev+1}, hmsh.deactivated{lev+1}))), elem);
    active{lev+1} = vertcat (active{lev+1}, new_possible_active_fun(new_functions));
    W{lev+1} = vertcat (W{lev+1}, zeros (sum (new_functions), 1));
  end
  
  if (isempty (W{lev}))
    active{lev} = zeros(0,1);
    W{lev} = zeros(0,1);
  end
end % for lev



% Update all the fields in hspace
if (~isempty(active{hspace.nlevels + 1}))
  hspace.nlevels = hspace.nlevels + 1;
end

hspace.active = active;
hspace.deactivated = deactivated;
hspace.coeff = cell2mat(W);
hspace.ndof_per_level = cellfun (@numel, hspace.active);
hspace.ndof = sum(hspace.ndof_per_level);

hspace.globnum_active = zeros (hspace.ndof,hmsh.ndim+1);
Nf = cumsum ([0; hspace.ndof_per_level(:)]);
for lev = 1:hspace.nlevels
  ind_f = (Nf(lev)+1):Nf(lev+1);
  if (~isempty(ind_f))
    hspace.globnum_active(ind_f, 1) = lev;
    globnum = cell (1,hmsh.ndim);
    [globnum{:}] = ind2sub (hspace.space_of_level(lev).ndof_dir, hspace.active{lev}(:));
    hspace.globnum_active(ind_f, 2:end) = cell2mat (globnum);
  end
end
