function [U, pts] = hspline_eval_old (u, hmsh, hspace, npts, varargin)
%
% function [U, pts] = hspline_eval (u, hmsh, hspace, npts, varargin)
%
% Evaluation of a hierarchical spline with degrees of freedom u
% 
% INPUT:    u: degrees of freedom
%           hmsh:
%           hspace:
%           npts: number of points for evaluation in each coordinate direction in each
%           cell. If npts == 0, the evaluation is performed in the
%           quadrature knots.
%           varargin:   'value', true or false, (Default: 'value', true)
%                       'gradient', true or false,
%                       'laplacian', true or false,
%
% OUTPUT:   U: values
%           pts: points
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ATENCION: Completar la descripcion de esta funcion
% 

if (hspace.ncomp ~= 1)
  disp('hspline_eval: Por ahora solo para escalares')
  return
end

if (nargin == 4)
  option = 'value';
elseif (nargin == 5)
  option = varargin{1};
end

C = hspace.C;

if (npts)
    npts_total = npts^hmsh.ndim;
else
    npts_total = hmsh.mesh_of_level(1).nqn;
end

if (strcmpi (option, 'value'))
    U = zeros (npts_total, hmsh.nel);
elseif (strcmpi (option, 'gradient'))
    U = zeros (hmsh.rdim, npts_total , hmsh.nel);
elseif (strcmpi (option, 'laplacian'))
    U = zeros (npts_total, hmsh.nel);
end

ndofs = 0;
nel = 0;
ndof_per_level = hspace.ndof_per_level;
dif = hmsh.nlevels - hspace.nlevels;
if dif
    ndof_per_level = [ndof_per_level(:); zeros(dif,1)];
end

pts = [];

for ilev = 1:hmsh.nlevels % Active levels
    nel_old = nel+1;
    nel = nel + hmsh.nel_per_level(ilev);
    ndofs = ndofs + ndof_per_level(ilev);
    dofs = 1:ndofs; % Active dofs from level 1 to ilev, in the numbering of the hierarchical space, whatever it is
    if (hmsh.msh_lev{ilev}.nel ~= 0)
        if npts
            rule = cell(hmsh.ndim,1);
            for idim = 1:hmsh.ndim
                rule{idim} = [linspace(-1+1e-12, 1-1e-12, npts); zeros(1, npts)];
            end
            qn = msh_set_quad_nodes (hmsh.mesh_of_level(ilev).breaks, rule);
            msh_plot = msh_cartesian (hmsh.mesh_of_level(ilev).breaks, qn, [], hmsh.geometry);
            sp_plot = hspace.space_of_level(ilev).constructor (msh_plot);
            
            msh_level = msh_evaluate_element_list (msh_plot, hmsh.active{ilev});
            sp_level = sp_evaluate_element_list (sp_plot, msh_level);
        else
            msh_level = hmsh.msh_lev{ilev};
            sp_level = hspace.sp_lev{ilev};
        end
        
        if nargout == 2
            pts = cat(3,pts, msh_level.geo_map);
        end
        
        if (strcmpi (option, 'value'))
            [eu, F] = sp_eval_msh_old (C{ilev}*u(dofs), sp_level, msh_level);
            U(:,nel_old:nel) = eu;
        elseif (strcmpi (option, 'gradient'))
            [eu, F] = sp_eval_grad_msh_old (C{ilev}*u(dofs), sp_level, msh_level);
            U(:,:,nel_old:nel) = eu;
        elseif (strcmpi (option, 'laplacian'))
            [eu, F] = sp_eval_lapl_msh_old (C{ilev}*u(dofs), sp_level, msh_level);
            U(:,nel_old:nel) = eu;
        end
        
    end
end
