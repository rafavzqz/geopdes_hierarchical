function fig = plot_numerical_and_exact_solution (u, hspace, geometry, npts, uex, fig)
%
% function plot_numerical_and_exact_solution(u, hspace, geometry, npts, uex)
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%
% Funcion personal
%
%
if (nargin == 6 && ishandle(fig))
  figure (fig)
elseif (nargout == 1)
  fig = figure;
else
  figure;
end

if (isa (hspace, 'hierarchical_space'))
  if (isa (hspace.space_of_level(1), 'sp_scalar'))
    sp_aux = hspace.space_of_level(1);
    is_vector = false;
  elseif (isa (hspace.space_of_level(1), 'sp_vector'))
    sp_aux = hspace.space_of_level(1).scalar_spaces{1};
    is_vector = true;
  end
elseif (isa (hspace, 'hierarchical_space_mp'))
  if (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_scalar'))
    sp_aux = hspace.space_of_level(1).sp_patch{1};
    is_vector = false;
  elseif (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_vector'))
    sp_aux = hspace.space_of_level(1).sp_patch{1}.scalar_spaces{1};
    is_vector = true;
  end
end
ndim = numel (sp_aux.knots);
first_knot = cellfun (@(x) x(1), sp_aux.knots);
last_knot = cellfun (@(x) x(end), sp_aux.knots);

if ((nargin > 4 && ~isempty (uex)))
  if (ndim == 3)
    warning ('The plot of the exact solution for volumes is not implemented. Plotting the numerical solution')
  elseif (is_vector)
    warning ('The plot of the exact solution for vectors is not implemented. Plotting the numerical solution')
  else
    subplot (1, 2, 2)
    for idim = 1:ndim
      x{idim} = linspace (first_knot(idim), last_knot(idim), npts(idim));
    end
    for iptc = 1:numel(geometry)
      F = geometry(iptc).map (x);
      rdim = size (F, 1);
      F = reshape (F, [rdim, npts]);

      if (ndim == 1 && rdim == 1)
        plot (F, uex(F)); hold on
      elseif (ndim == 1 && rdim == 2)
        plot3 (F(1,:), F(2,:), uex(F)); hold on
      elseif (ndim == 2 && rdim == 2)
        surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(uex(F(1,:,:),F(2,:,:)))); hold on
        shading interp
        xlabel('x'); ylabel('y'); zlabel('z')
      elseif (ndim == 2 && rdim == 3)
        surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)), squeeze(uex(F(1,:,:),F(2,:,:),F(3,:,:)))); hold on
        shading interp
        xlabel('x'); ylabel('y'); zlabel('z')
      else
        error('The plot of the exact solution with ndim=%d and rdim=%d is not implemented', ndim, rdim)
      end
    end
    title ('Exact solution')
    hold off
    subplot (1, 2, 1)
  end
end
sp_plot_solution (u, hspace, geometry, npts)
shading interp
title ('Numerical solution'),

end
