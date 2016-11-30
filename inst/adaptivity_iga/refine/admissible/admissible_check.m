function [flag] = admissible_check(hmsh, hspace, m)
%returns 1 if the hierarchical mesh is admissible of class m, 0 otherwise

support_active=cell(hmsh.nlevels,hspace.ndof);
active_cell=cell(1,hmsh.nel);

%lev=1;
ndof_prev_levs=0;
for j=1:hspace.ndof
%     while hspace.ndof_per_level(lev)==0 && lev<hspace.nlevels
%         lev=lev+1;
%     end
%     jl=hspace.active{lev}(j-ndof_prev_levs);
    for  le=hmsh.nlevels
        support_active(:,j)=supportTHB(hmsh, hspace, j);
    end
%     while j>=ndof_prev_levs+hspace.ndof_per_level(lev) && lev<hspace.nlevels
%         ndof_prev_levs=ndof_prev_levs+hspace.ndof_per_level(lev);
%         lev=lev+1;
%     end
end

lev=1;
nel_prev_levs=0;
for i=1:hmsh.nel
    while hmsh.nel_per_level(lev)==0 && lev<hmsh.nlevels
        lev=lev+1;
    end
    il=hmsh.active{lev}(i-nel_prev_levs);
    levf=1;
    ndof_prev_levs=0;
    for j=1:hspace.ndof
        while hspace.ndof_per_level(levf)==0 && levf<hspace.nlevels
            levf=levf+1;
        end
        if isempty(find(support_active{lev,j}==il,1))==0
            active_cell{i}=[active_cell{i} levf];
        end
        while j>=ndof_prev_levs+hspace.ndof_per_level(levf) && levf<hspace.nlevels
            ndof_prev_levs=ndof_prev_levs+hspace.ndof_per_level(levf);
            levf=levf+1;
        end
    end 
    active_cell{i}=unique(active_cell{i}); 
    while i>=nel_prev_levs+hmsh.nel_per_level(lev) && lev<hmsh.nlevels
        nel_prev_levs=nel_prev_levs+hmsh.nel_per_level(lev);
        lev=lev+1;
    end
end

flag=1;
i=1;
while (flag==1) && (i<=hmsh.nel)
    if length(active_cell{i})>=m+1
        flag=0;
    else
        i=i+1;
    end
end

end