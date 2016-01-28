function plot_numerical_and_exact_solution(u, hspace, geometry, npts, uex)
%
% function plot_numerical_and_exact_solution(u, hspace, geometry, npts, uex)
% 
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% 
% Funcion personal
% 
%

ndim = numel (npts);
if (nargin < 5)
  uex = [];
end

switch ndim
  case 1, %x = hmsh.quad_nodes(:);
    [Z, pts] = sp_eval (u, hspace, geometry, npts);
    if (isempty (uex))
      figure
    else
      subplot (1,2,2)
      plot (pts, uex(pts));
      title('Exact solution'),
      subplot (1,2,1)
    end
    plot (pts, Z)
    title ('Numerical solution'),
        
  case 2,
    [Z, pts] = sp_eval (u, hspace, geometry, npts);
    x = reshape (pts(1,:,:), npts);
    y = reshape (pts(2,:,:), npts);
    if (isempty (uex))
      figure
    else
      subplot (1,2,2)
      surf (x, y, uex(x,y));
      shading interp
      title('Exact solution'),
      axis tight
      xlabel('x')
      ylabel('y')
      hold off
      subplot (1,2,1)
    end
    surf (x, y, Z);
    shading interp
    title ('Numerical solution'),
    axis tight
    xlabel('x')
    ylabel('y')
    hold off
        
    case 3, disp('Not implemented'),
end
