function [list_of_cells] = support_ext (hmsh, hspace, Q_ind, lev_Q, lev_s)
%SUPPORT_EXT returns the list of indices of the cells of tensor-product mesh
%of level lev_s which belong to the support extension of the cell(s) of index 
%Q_ind and level lev_Q

if (lev_Q == lev_s)
    funs = sp_get_basis_functions (hspace.space_of_level(lev_Q), hmsh.mesh_of_level(lev_Q), Q_ind);
    list_of_cells = sp_get_cells (hspace.space_of_level(lev_Q), hmsh.mesh_of_level(lev_Q), funs);
else
    ancestors = get_ancestors (hmsh, Q_ind, lev_Q, lev_s);
    list_of_cells = support_ext (hmsh, hspace, ancestors, lev_s, lev_s);
end

end