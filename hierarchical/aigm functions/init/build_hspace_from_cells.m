function [hmsh, hspace] = build_hspace_from_cells(dim, p, initial_num_el, cells,space_type, boundary, graficar_malla)
%
% function [hmsh, hspace] = build_hspace_from_cells(dim, p, initial_num_el, cells, space_type, boundary, graficar_malla)
%
% This function fills hmsh and hspace. The active cells are given in
% cells{lev}, for lev = 1,2,...
%
% ATENCION: Completar la descripcion de esta funcion y decidir si queremos argumentos de entrada mas generales. Por ahora, es una construccion
% particular
%

if nargin == 6
    graficar_malla = 1;
end

switch dim
    case 1, problem_data.geo_name = nrbline ([0 0], [1 0]);
    case 2, problem_data.geo_name = 'geo_square.txt'; 
    case 3, problem_data.geo_name = 'geo_cube.txt';
end

% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree     = p*ones(1,dim);       % Degree of the splines
method_data.regularity = method_data.degree-1;       % Regularity of the splines
method_data.nsub       = initial_num_el*ones(1,dim);       % Number of subdivisions
method_data.nquad      = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.space_type = space_type;           % 0: , 1: Full basis (B-splines)

[hmsh, hspace] = adaptivity_initialize (problem_data, method_data);

nref = numel(cells);

for ref = 1:nref
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MARK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    marked = cell(ref,1);
    marked(ref) = cells(ref);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REFINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hmsh, hspace] = refine (hmsh, hspace, marked, 'elements', boundary);
    
end

if graficar_malla
plot_msh(hmsh)
end
