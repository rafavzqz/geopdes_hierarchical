function [list_of_cells] = supportTHB(hmsh, hspace, B_ind)
%returns a cell array containing, for each level of the hierarchical mesh,
%the list of active cells belonging to the support of the THB-spline

lev=1;
ndof_prev_levs=0;
while B_ind>ndof_prev_levs+hspace.ndof_per_level(lev) && lev<hspace.nlevels
    ndof_prev_levs=ndof_prev_levs+hspace.ndof_per_level(lev);
    lev=lev+1;
end
lev_B=lev;
%at this point we computed the level of the THB-spline with (global) index B_ind

list_of_functions=cell(hspace.nlevels,1);
list_of_cells=cell(hmsh.nlevels,1);

%initialize the cell-array containing, for each level, the B-spline to be
%considered composing the support of the THB-spline
for l=1:lev_B
    list_of_functions{l}=[];
end
for l=lev_B+1:hspace.nlevels
    list_of_functions{l}=setdiff(find(hspace.Csub{l}(:,B_ind)),union(hspace.active{l},hspace.deactivated{l}));
end

for l=1:hspace.nlevels
    list_of_cells{l}=[];
end

for h=lev_B+1:hspace.nlevels
    [aux,~]=sp_get_cells(hspace.space_of_level(h), hmsh.mesh_of_level(h), list_of_functions{h});
    for l=lev_B:hspace.nlevels
        list_of_cells{l}=[list_of_cells{l}; ancestry_new(hmsh,aux,h,l)];
        list_of_cells{l}=[list_of_cells{l}; progeny_new(hmsh,aux,h,l)];
        list_of_cells{l}=unique(list_of_cells{l});
        list_of_cells{l}=intersect(list_of_cells{l},hmsh.active{l});
    end
end

end