function [ list_of_cells ] = neigh(hmsh, hspace, Q_ind, lev_Q, m)
%NEIGH computes the neighborhood of the cell(s) Q_ind of level lev_Q belonging to
%the hierarchical mesh hmsh (admissible of class m) associated to the space hspace,
%according to definition in of [A. Buffa, C. Giannelli, 
%Adaptive isogeometric methods with hierarchical splines: Error estimator and convergence, 2016]

list_of_cells = [];
lev_s = lev_Q-m+2;
if lev_s <= 1
    list_of_cells = [];
else
    all_el_lev = 1:hmsh.mesh_of_level(lev_s).nel;
    inactive_el = setdiff (all_el_lev,hmsh.active{lev_s});
    el = intersect (inactive_el,support_ext(hmsh,hspace,Q_ind,lev_Q, lev_s));
    anc = get_ancestors (hmsh,el,lev_s,lev_s-1);
    list_of_cells = union (list_of_cells, intersect(anc,hmsh.active{lev_s-1}));
end

end