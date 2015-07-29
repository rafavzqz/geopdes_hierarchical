function [errh1, errl2, errh1s, errh1_elem, errl2_elem, errh1s_elem] = compute_H1_error(u, uex, graduex, hmsh, hspace)
%
% function [errh1, errl2, errh1s, errh1_elem, errl2_elem, errh1s_elem] = compute_H1_error(u, uex, graduex, hmsh, hspace)
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ATENCION: Mejorar y documentar esta funcion
%

if hspace.ncomp ~= 1
    disp('compute_H1_error: Por ahora solo para escalares')
    return
end

[val, pts] = hspline_eval (u, hmsh, hspace, 0);
valgrad = hspline_eval (u, hmsh, hspace, 0, 'gradient');

switch hmsh.ndim
    case 1,
        valex = uex (pts(1,:,:));
        gradex = graduex (pts(1,:,:));
    case 2,
        valex = uex (pts(1,:,:), pts(2,:,:));
        gradex = graduex (pts(1,:,:), pts(2,:,:));
    case 3,
        valex = uex (pts(1,:,:), pts(2,:,:), pts(3,:,:));
        gradex = graduex (pts(1,:,:), pts(2,:,:), pts(3,:,:));
end

nqn = hmsh.mesh_of_level(1).nqn;
valex = reshape (valex, nqn, hmsh.nel);
gradex = reshape (gradex, hmsh.rdim, nqn, hmsh.nel);

% jw = hmsh.quad_weights .* hmsh.jacdet;
quad_weights = [];
jacdet = [];
for ilev = 1:hmsh.nlevels % Active levels
    if (hmsh.msh_lev{ilev}.nel ~= 0)
            quad_weights = cat(2,quad_weights, hmsh.msh_lev{ilev}.quad_weights);
            jacdet = cat(2,jacdet, hmsh.msh_lev{ilev}.jacdet);
    end
end
jw = quad_weights .* jacdet;


errl2_elem = sqrt (sum (jw .* (val - valex).^2, 1));
errl2 = norm (errl2_elem);

errh1s_elem = sqrt (sum (reshape (sum ((valgrad - gradex).^2, 1), nqn, hmsh.nel) .*jw, 1));

errh1s = norm(errh1s_elem);

errh1_elem = sqrt (errl2_elem.^2 + errh1s_elem.^2);
errh1 = norm (errh1_elem);