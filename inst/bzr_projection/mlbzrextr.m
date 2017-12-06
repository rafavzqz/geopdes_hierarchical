function [C_mle, hier_funct, Cext_e, C_ml] = mlbzrextr(hmsh, hspace, active, lev, index)

space = hspace.space_of_level(lev);
msh = hmsh.mesh_of_level(lev);

knots = space.knots;
degree = space.degree;
funct = sp_get_basis_functions(space, msh, index);
[I,J,K] = ind2sub(msh.nel_dir, index);
sub = [I; J; K];
Cext_e = 1;
for iDir = 1:hmsh.ndim
    C_dir = bzrextr(knots{iDir},degree(iDir));
    Cext_e = kron(C_dir(:,:,sub(iDir)), Cext_e);
end

ndof_per_level = cellfun (@numel, active);
hier_funct = [];

for iLev = 1:lev-1
    f = intersect(active{iLev}, sp_get_basis_functions(hspace.space_of_level(iLev), hmsh.mesh_of_level(iLev), get_ancestors(hmsh,index,lev,iLev)));
    hier_funct = [hier_funct; find(ismember(active{iLev},f)) + sum(ndof_per_level(1:iLev-1))];
end


f = intersect(active{lev}, funct);
hier_funct = [hier_funct; find(ismember(active{lev},f)) + sum(ndof_per_level(1:lev-1))];

C_ml = transpose(hspace.Csub{lev}(funct, hier_funct));
C_mle = C_ml * Cext_e;

end