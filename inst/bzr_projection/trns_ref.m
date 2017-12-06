function [L] = trns_ref(hmsh, hspace, active_trg, deactivated_trg, lev, index)

space = hspace.space_of_level(1);
msh = hmsh.mesh_of_level(1);

if lev > 1
    base_msh_el = get_ancestors(hmsh,index,lev,1);
    fct_base = sp_get_basis_functions(space,msh,base_msh_el);
    deactivated = intersect(deactivated_trg{1}, fct_base);
    active_0 = find(ismember(fct_base, active_trg{1}));
    deactivated_0 = find(ismember(fct_base, deactivated_trg{1}));
    n_active_0 = length(active_0);
    n_deactivated = length(deactivated_0);
    n_active_deactivated_0 = n_active_0 + n_deactivated;
    L = eye(length(fct_base),length(fct_base));
    J = L(deactivated_0,deactivated_0);
    L(deactivated_0, :) = zeros(n_deactivated, length(fct_base));
    L(:, deactivated_0) = [];
    n_active = n_active_0;
    ext_active = active_0;
else
    fct_base = sp_get_basis_functions(space,msh,index);
    L = eye(length(fct_base),length(fct_base));
end

for iLev = 2:lev-1
    space = hspace.space_of_level(iLev);
    msh = hmsh.mesh_of_level(iLev);
    parent = get_ancestors(hmsh,index,lev,iLev);
    fct_support = sp_get_basis_functions(space,msh,parent);
    active = intersect(active_trg{iLev}, fct_support);
    
    aux = matrix_basis_change__ (hspace, iLev);
    J = J*aux(:,deactivated)';
    L_temp = zeros(n_active_deactivated_0, n_active + length(active));
    L_temp(1:size(L,1), 1:numel(ext_active)) = L;
    if ~isempty(active)
        L_temp(deactivated_0, n_active+1:end) = J(:, active);
    end
    L = L_temp;
    deactivated = intersect(deactivated_trg{iLev}, fct_support);
    n_active = n_active + length(active);
    J = J(:,deactivated);
    ext_active = [ext_active; active];
end

if lev > 1
    space = hspace.space_of_level(lev);
    msh = hmsh.mesh_of_level(lev);
    fct_support = sp_get_basis_functions(space,msh,index);
    active = intersect(active_trg{lev}, fct_support);
    aux =  matrix_basis_change__ (hspace, lev);
    J = J*aux(:,deactivated)';
    L_temp = zeros(n_active_deactivated_0, n_active + length(active));
    L_temp(1:size(L,1), 1:numel(ext_active)) = L;
    L_temp(deactivated_0, n_active+1:end) = J(:, active);
    L = L_temp;
end


end


