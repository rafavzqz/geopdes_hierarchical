function est = estimate (u, hmsh, hspace, problem_data, est_type)

% INPUT:
%
%     hspace:
%     hmsh:
%     u:        vector of dof weights
%     problem_data.f:     function handle for the rhs function

if strcmp(est_type,'none')
    est = ones(1,hspace.ndof)/sqrt(hspace.ndof);
    return;
end

% pts = [];
% for ilev = 1:hmsh.nlevels % Active levels
%     if (hmsh.msh_lev{ilev}.nel ~= 0)
%             pts = cat(3,pts, hmsh.msh_lev{ilev}.geo_map);
%     end
% end

[der2num, pts] = hspline_eval(u,hmsh,hspace,0,'laplacian'); %  ones( hmsh.nqn, hmsh.nel);

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

% w = hmsh.quad_weights .* hmsh.jacdet;
quad_weights = [];
jacdet = [];
for ilev = 1:hmsh.nlevels % Active levels
    if (hmsh.msh_lev{ilev}.nel ~= 0)
            quad_weights = cat(2,quad_weights, hmsh.msh_lev{ilev}.quad_weights);
            jacdet = cat(2,jacdet, hmsh.msh_lev{ilev}.jacdet);
    end
end
w = quad_weights .* jacdet;


[h, ms] = get_meshsize(hmsh);

switch est_type
    case 'estandar',
        est = sqrt (sum (aux.*w));
        est = h.*est(:);
    case 'nonestandar',
        coef = ms(hspace.globnum_active(:,1)).^2.*hspace.coeff(:);
        H2 = hspline_eval(coef, hmsh, hspace,0);
        est = aux.*H2;
        est = sqrt (sum (est.*w));
    case 'basis_functions',
        forma = 1; % 0 o 1
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
        
%         est1 =  zeros(hspace.ndof,1);
%         for i = 1:hspace.ndof
%             uu = zeros(hspace.ndof,1);
%             uu(i) = 1;
%             basis_fun_val = hspline_eval(uu, hmsh, hspace);
%             if forma
%                 est1(i) =coef(i)*sqrt(sum(sum((aux.*basis_fun_val).*w)));
%             else
%                 lev = hspace.globnum_active(i,1);
%                 cells = get_cells(hspace.globnum_active(i,2:end), hspace.degree, hmsh.mesh_of_level(lev).nel_dir);
%                 supp_size = max(cells)-min(cells)+1;
%                 diametro = norm(supp_size./hmsh.mesh_of_level(lev).nel_dir);
%                 est1(i)  = diametro*coef(i)*sqrt(sum(sum((aux.*basis_fun_val).*w)))/pi;
%             end
%         end
    case 'basis_functions_2'
        [d,diametro] = distance_to_boundary(hmsh, hspace, squeeze(pts(1,:,:)), squeeze(pts(2,:,:)));
        est = zeros(hspace.ndof,1);
        for i = 1:hspace.ndof
            uu = zeros(hspace.ndof,1);
            uu(i) = 1;
            basis_fun_val = hspline_eval(uu, hmsh, hspace,0);
            est(i)  = diametro(i)*sqrt(hspace.coeff(i)*sum(sum((aux.*d{i}.^2.*basis_fun_val).*w)));
        end
        
end
