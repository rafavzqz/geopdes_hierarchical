function est = estimate_laplace_param (u, hmsh, hspace, problem_data, adaptivity_data)
%
% function est = estimate_laplace_param (u, hmsh, hspace, problem_data, adaptivity_data)
%
% This function computes some local a posteriori error estimators for the
% laplacian (problem_data.c_diff = 1), where the vector u contains the degrees of freedom of the Galerkin solution
%
% INPUT:
%
%     hspace:
%     hmsh:
%     u:        degrees of freedom
%     problem_data.f:     function handle for the rhs function
%     adaptivity_data.flag: 'elements' or 'functions' or 'none'
%
% OUTPUT:  est: array with the values of the local error estimators
%
% This function uses:   get_meshsize_param
%                       hspline_eval
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ATENCION: Deberiamos modificar get_meshsize_param para poder calcular los
% estimadores en un dominio fisico distinto del parametrico. 
% Solo hay que calcular el meshsize del dominio fisico, porque las
% evaluaciones ya se estan haciendo en los nodos de cuadratura asi que eso
% no sera necesario cambiarlo.
%

if strcmp(adaptivity_data.flag,'none')
    est = ones(1,hspace.ndof)/sqrt(hspace.ndof);
    return;
end

[der2num, pts] = hspline_eval(u,hmsh,hspace,0,'laplacian'); 

switch hmsh.ndim
    case 1,
        valf = feval (problem_data.f, pts);
    case 2,
        valf = problem_data.f (pts(1,:,:), pts(2,:,:));
    case 3,
        valf = problem_data.f (pts(1,:,:), pts(2,:,:), pts(3,:,:));
end

valf = squeeze(valf);

aux = (valf + der2num).^2; % size(aux) = [hmsh.nqn, hmsh.nel], valores en los nodos de cuadratura

quad_weights = [];
jacdet = [];
for ilev = 1:hmsh.nlevels % Active levels
    if (hmsh.msh_lev{ilev}.nel ~= 0)
            quad_weights = cat(2,quad_weights, hmsh.msh_lev{ilev}.quad_weights);
            jacdet = cat(2,jacdet, hmsh.msh_lev{ilev}.jacdet);
    end
end
w = quad_weights .* jacdet;

[h, ms] = get_meshsize_param(hmsh);

switch adaptivity_data.flag
    case 'elements',
        est = sqrt (sum (aux.*w));
        est = h.*est(:);
    case 'functions',
        forma = 1; % 0 o 1
        % La version final sera con la opcion forma = 0 directamente. Ahora
        % estoy realizando comparaciones con corridas anteriores en las que
        % forma era igual a 1.
        if forma
            coef = ms(hspace.globnum_active(:,1)).*sqrt(hspace.coeff(:));
        else
            coef = sqrt(hspace.coeff(:));
        end
        
        est = zeros(hspace.ndof,1);
        ndof_per_level = hspace.ndof_per_level;
        dif = hmsh.nlevels - hspace.nlevels;
        if dif
            ndof_per_level = [ndof_per_level(:); zeros(dif,1)];
        end
        ndofs = 0;
        Ne = cumsum([0; hmsh.nel_per_level(:)]);
        for ilev = 1:hmsh.nlevels % Active levels
            ndofs = ndofs + ndof_per_level(ilev);
            if hmsh.msh_lev{ilev}.nel
                ind_e = (Ne(ilev)+1):Ne(ilev+1);
                b_lev = op_f_v (hspace.sp_lev{ilev}, hmsh.msh_lev{ilev}, aux(:,ind_e));
                
                dofs = 1:ndofs; % Active dofs from level 1 to ilev, in the numbering of the hierarchical space, whatever it is
                est(dofs) = est(dofs) + hspace.C{ilev}'*b_lev;
            end
        end
        if forma
            est = coef.*sqrt(est);
        else
            diametro = zeros(hspace.ndof,1);
            for i = 1:hspace.ndof
                lev = hspace.globnum_active(i,1);
                cells = get_cells(hspace.globnum_active(i,2:end), hspace.degree, hmsh.mesh_of_level(lev).nel_dir);
                supp_size = max(cells)-min(cells)+1;
                diametro(i) = norm(supp_size./hmsh.mesh_of_level(lev).nel_dir);
            end
            est = diametro.*coef.*sqrt(est)/pi;
        end
       
        
end
