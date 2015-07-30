% HMSH_PLOT_CELLS: plot the cells of the hierarchical mesh
%
%   hmsh_plot_cells (hmsh, [fig_number])
%
% INPUT:
%
%    hmsh:       object representing the hierarchical mesh (see hierarchical_mesh)
%    fig_number: figure number where to plot (if not given, a new figure is open)
%
% The function is still unefficient, but it can be useful
%
% Copyright (C) 2015 Eduardo M. Garau, Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function hmsh_plot_cells (hmsh, nfig)

if (nargin == 1)
  figure
  hold_flag = 0;
else
  figure(nfig)
  hold_flag = ishold;
end
hold on

npts = 2;
for ilev = 1:hmsh.nlevels
  if (hmsh.nel_per_level(ilev) > 0)
    rule = cell(hmsh.ndim,1);
    for idim = 1:hmsh.ndim
      rule{idim} = [linspace(-1+1e-12, 1-1e-12, npts); zeros(1, npts)];
    end
    qn = msh_set_quad_nodes (hmsh.mesh_of_level(ilev).breaks, rule);
    msh_plot = msh_cartesian (hmsh.mesh_of_level(ilev).breaks, qn, [], hmsh.geometry);
    msh_level = msh_evaluate_element_list (msh_plot, hmsh.active{ilev});

    x = cell (hmsh.rdim, 1);
    for idim = 1:hmsh.rdim
      x{idim} = reshape (msh_level.geo_map(idim,:,:), [npts * ones(1, hmsh.ndim), msh_level.nel]);
    end
    for idim = hmsh.rdim+1:3
      x{idim} = zeros (size (x{1}));
    end
    
    if (hmsh.ndim == 1)
      plot3 (x{1}, x{2}, x{3}, 'k', 'Marker', 'x');
    elseif (hmsh.ndim == 2)
      for iel = 1:msh_level.nel
        surf (x{1}(:,:,iel), x{2}(:,:,iel), x{3}(:,:,iel));
      end
    elseif (hmsh.ndim == 3)
      for iel = 1:msh_level.nel
        siz = [2 2];
        surf (reshape (x{1}(1,:,:,iel), siz), reshape (x{2}(1,:,:,iel), siz), reshape (x{3}(1,:,:,iel), siz), 'FaceAlpha', 0);
        surf (reshape (x{1}(:,1,:,iel), siz), reshape (x{2}(:,1,:,iel), siz), reshape (x{3}(:,1,:,iel), siz), 'FaceAlpha', 0);
        surf (reshape (x{1}(:,:,1,iel), siz), reshape (x{2}(:,:,1,iel), siz), reshape (x{3}(:,:,1,iel), siz), 'FaceAlpha', 0);
        surf (reshape (x{1}(end,:,:,iel), siz), reshape (x{2}(end,:,:,iel), siz), reshape (x{3}(end,:,:,iel), siz), 'FaceAlpha', 0);
        surf (reshape (x{1}(:,end,:,iel), siz), reshape (x{2}(:,end,:,iel), siz), reshape (x{3}(:,end,:,iel), siz), 'FaceAlpha', 0);
        surf (reshape (x{1}(:,:,end,iel), siz), reshape (x{2}(:,:,end,iel), siz), reshape (x{3}(:,:,end,iel), siz), 'FaceAlpha', 0);
      end
    end
  end
end


if (~hold_flag)
  hold off
end

end
