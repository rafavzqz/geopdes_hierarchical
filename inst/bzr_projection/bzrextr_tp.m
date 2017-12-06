function [C] = bzrextr_tp(hmsh, hspace, lev)

space = hspace.space_of_level(lev);
msh = hmsh.mesh_of_level(lev);

knots = space.knots;
degree = space.degree;
C = ones(prod(degree+1),prod(degree+1),hmsh.nel);

for el=1:msh.nel
    C_el = 1;
    [I,J,K] = ind2sub(msh.nel_dir, el);
    sub = [I; J; K];
    for iDir = 1:hmsh.ndim
        p = degree(iDir);
        C_e = bzrextr(knots{iDir},p);
        C_el = kron(C_e(:,:,sub(iDir)),C_el);
        
    end
    C(:,:,el) = C_el;
end

end

