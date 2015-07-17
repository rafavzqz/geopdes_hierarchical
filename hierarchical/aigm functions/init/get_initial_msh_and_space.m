function [msh, space] = get_initial_msh_and_space(dim, p, initial_num_el, geometry)
%
% function [msh, space] = get_initial_msh_and_space(dim, p, initial_num_el, geometry)
%
% Completar descripcion
%

% knots = geometry.nurbs.knots;
degree = p*ones(1,dim);
%degree = [2 3];
nelem =  initial_num_el*ones(1,dim);
break_points = cell(1,dim);
open_knot_vector = cell(1,dim);
nbreaks = nelem+1;
for idir = 1:dim
    break_points{idir} = linspace(0,1,nbreaks(idir));
    open_knot_vector{idir} = [ zeros(1, degree(idir)) linspace(0,1,nbreaks(idir)) ones(1, degree(idir))];
end
[qn, qw] = msh_set_quad_nodes (break_points, msh_gauss_nodes(degree+1));
msh = msh_cartesian (break_points, qn, qw, geometry);
space = sp_bspline(open_knot_vector, degree, msh);