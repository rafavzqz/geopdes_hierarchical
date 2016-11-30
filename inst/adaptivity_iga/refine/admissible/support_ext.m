function [list_of_cells] = support_ext(hmsh, hspace, Q_ind, lev_Q, lev_s)
%SUPPORT_EXT returns the list of indices of the cells of tensor-product mesh
%of level lev_s which belong to the support extension of the cell(s) of index 
%Q_ind and level lev_Q

dd=hspace.space_of_level(lev_Q).degree; %degree vector
n=length(dd); %dimension

if lev_Q==lev_s
    subindices = cell (n, 1);
    [subindices{:}]=ind2sub ([hmsh.mesh_of_level(lev_Q).nel_dir, 1], Q_ind);
    %getting indices (in each direction) of the cells Q_ind
    list_of_cells_per_cell = cell (numel(Q_ind), 1);
    for icell = 1:numel(Q_ind)
        Qa=cell(n,1);
        Qb=cell(n,1);
        for h=1:n
            %h-th coordinate of the "left-bottom" and "right-top" vertices of the cell icell
            Q_a{h}=[hmsh.mesh_of_level(lev_Q).breaks{h}(subindices{h}(icell))]; 
            Q_b{h}=[hmsh.mesh_of_level(lev_Q).breaks{h}(subindices{h}(icell)+1)]; 
            %h-th coordinate of the "left-bottom" and "right-top" vertices of the support extension of the cell icell...
            Q_knota_ind{h}=find(hspace.space_of_level(lev_Q).knots{h}==Q_b{h},1)-dd(h)-1;
            Q_knotb_ind{h}=find(hspace.space_of_level(lev_Q).knots{h}==Q_a{h},1,'last')+dd(h)+1;
            %...and corresponding indices in the knot vector
            Q_knota{h}=hspace.space_of_level(lev_Q).knots{h}(Q_knota_ind{h});
            Q_knotb{h}=hspace.space_of_level(lev_Q).knots{h}(Q_knotb_ind{h});
            %h-th index of the "left-bottom" and "right-top" elements of the support extension
            Q_as{h}=find(hmsh.mesh_of_level(lev_Q).breaks{h}==Q_knota{h});
            Q_bs{h}=find(hmsh.mesh_of_level(lev_Q).breaks{h}==Q_knotb{h})-1;
            cells_1d{h}=[Q_as{h}:Q_bs{h}];
        end
        cells=cell (n, 1);
        [cells{:}] = ndgrid (cells_1d{:});
        %list of cells  of the support extension of the cell icell
        list_of_cells_per_cell{icell}= sub2ind ([hmsh.mesh_of_level(lev_Q).nel_dir, 1], cells{:});
    end
    list_of_cells_per_cell = cellfun(@(x) x(:), list_of_cells_per_cell, 'UniformOutput', false);
    %list of cells  of the support extension of all the cells cells Q_ind
    list_of_cells = unique (vertcat (list_of_cells_per_cell{:}));
else
    list_of_cells=support_ext(hmsh,hspace,get_ancestors(hmsh,Q_ind,lev_Q,lev_s),lev_s,lev_s);
end

end