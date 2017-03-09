%% EX_ADAPTIVITY_THBCOARSENING2D
% Test for THB-spline coarsening algorithm:
% A 3 level truncated hierarchical 2D mesh is considered, one level coarsening
% is performed sequentially between the second and the first level and between the third
% and the second level. N.B. The method works ONLY if THB-splines are used.
%
% Copyright (C) 2017 Massimo Carraturo
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

clear problem_data
close all
clc
% Initial domain, defined as NURBS map given in a text file
problem_data.geo_name = 'geo_square.txt';

% CHOICE OF THE DISCRETIZATION PARAMETERS (Initial mesh)
clear method_data
method_data.degree      = [2 2];        % Degree of the splines
method_data.regularity  = [1 1];        % Regularity of the splines
method_data.nsub_coarse = [2 2];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad       = [3 3];        % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';   % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = 1;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';

% GRAPHICS
plot_data.npoints_x = 100;       %number of points x-direction in post-processing
plot_data.npoints_y = 100;       %number of points y-direction in post-processing
plot_data.npoints_z = 1;         %number of points z-direction in post-processing

%% Generate the first level of initial mesh
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (problem_data.geo_name);
[knots, zeta] = kntrefine (geometry.nurbs.knots, method_data.nsub_coarse-1, method_data.degree, method_data.regularity);
% tensor product level 1
rule     = msh_gauss_nodes (method_data.nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh   = msh_cartesian (zeta, qn, qw, geometry);
space = sp_bspline (knots, method_data.degree, msh);
% hierarchical space
hmsh     = hierarchical_mesh (msh, method_data.nsub_refine);
hspace   = hierarchical_space (hmsh, space, method_data.space_type, method_data.truncated);

%% Add second and thrid level using THB-refinement
% Second level
marked_ref = cell(1, hspace.nlevels);
marked_ref{1} = [1 2 3 4];
[hmsh, hspace, ~] = adaptivity_refine (hmsh, hspace, marked_ref, adaptivity_data);
hmsh_plot_cells (hmsh, 20, 1 );

% add new level

hmsh.nlevels = hmsh.nlevels + 1;
hmsh.active{hmsh.nlevels} = [];
hmsh.deactivated{hmsh.nlevels} = [];
hmsh.nel_per_level(hmsh.nlevels) = 0;
hmsh.mesh_of_level(hmsh.nlevels) = msh_refine (hmsh.mesh_of_level(hmsh.nlevels-1), hmsh.nsub);
hmsh.msh_lev{hmsh.nlevels} = [];


% assign dofs
v1 = [1 1.0 1.0 1.0 1 1 1 1 1 1]';
v2 = [1.0 1.1 1.1 1.1 1.3 1.4]';
v3 = [1.75 1.4 0.8 1.5 1.75 1.6]';
v4 = v3*1.5;

hspace.dofs = cat(1, v1, v2, v3, v4, v2, [1.4, 2]');
initial_values = hspace.dofs;

% plot initial state
npts = [plot_data.npoints_x plot_data.npoints_y ];
[eu, F] = sp_eval (initial_values, hspace, geometry, npts);
figure(11); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)

%% Refinement
marked_ref{1} = [];
marked_ref{2} = [6 7 10 11];
marked_ref{3} = [];
[hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked_ref, adaptivity_data);
hmsh_plot_cells (hmsh, 20, (figure(3)));
u_ref = Cref*initial_values;

marked_ref{1} = [];
marked_ref{2} = [];
marked_ref{3} = [28 29 36 37];
marked_ref{4} = [];
[hmsh, hspace, Cref] = adaptivity_refine (hmsh, hspace, marked_ref, adaptivity_data);
hmsh_plot_cells (hmsh, 20, (figure(3)));
hspace.dofs = Cref*u_ref;

% plot refined state
npts = [plot_data.npoints_x plot_data.npoints_y ];
[eu, F] = sp_eval (hspace.dofs, hspace, geometry, npts);
figure(12); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)

%% Coarsening back to the initial state

marked_coarse{1} = [];
marked_coarse{2} = [];
marked_coarse{3} = [];
marked_coarse{4} = cat(2,[103 104 119 120],[103+32 104+32 119+32 120+32],...
    [103+2 104+2 119+2 120+2], [103+34 104+34 119+34 120+34]);
[hmsh, hspace, u] = adaptivity_coarsen(hmsh, hspace, marked_coarse, adaptivity_data);
hspace.dofs = u;

marked_coarse{1} = [];
marked_coarse{2} = [];
marked_coarse{3} = [19 20 27 28 19+16 20+16 27+16 28+16];
[hmsh, hspace, u] = adaptivity_coarsen(hmsh, hspace, marked_coarse, adaptivity_data);
hmsh_plot_cells (hmsh, 20, (figure(4)));

% plot coarse state
npts = [plot_data.npoints_x plot_data.npoints_y];
[eu, F] = sp_eval (u, hspace, geometry, npts);
figure(13); surf (squeeze(F(1,:,:)), squeeze(F(2,:,:)), eu)

%% Check

for i=1:numel(u)
    if (u(i) < initial_values(i)-1.0e-08|| u(i) > initial_values(i)+1.0e-08)
        disp('is not a projector !');
        break;
    else
        disp('O.K.');
    end
end

