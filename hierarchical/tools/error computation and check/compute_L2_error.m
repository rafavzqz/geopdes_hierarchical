function [err, err_elem] = compute_L2_error(u, uex, hmsh, hspace)
%
% function [err, err_elem] = compute_L2_error(u, uex, hmsh, hspace)
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ATENCION: Mejorar y documentar esta funcion
%

if hmsh.ndim ~= 2
    disp('compute_L2_error: Por ahora solo para 2d')
    return
end

[val,pts] = hspline_eval(u, hmsh, hspace, 0);
valex = uex(pts(1,:,:), pts(2,:,:));
valex = reshape (valex, hmsh.nqn, hmsh.nel);
%jw = hmsh.quad_weights.*hmsh.jacdet;
quad_weights = [];
jacdet = [];
for ilev = 1:hmsh.nlevels % Active levels
    if (hmsh.msh_lev{ilev}.nel ~= 0)
            quad_weights = cat(2,quad_weights, hmsh.msh_lev{ilev}.quad_weights);
            jacdet = cat(2,jacdet, hmsh.msh_lev{ilev}.jacdet);
    end
end
jw = quad_weights .* jacdet;


err_elem = sqrt(sum (jw .* (val - valex).^2));
err = norm(err_elem);