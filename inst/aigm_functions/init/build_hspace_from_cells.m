function [hmsh, hspace] = build_hspace_from_cells(dim, p, initial_num_el, cells,space_type, graficar_malla, truncated)
%
% function [hmsh, hspace] = build_hspace_from_cells(dim, p, initial_num_el, cells, space_type, graficar_malla, truncated)
%
% This function fills hmsh and hspace. The active cells are given in
% cells{lev}, for lev = 1,2,...
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%
% ATENCION: 
% - Decidir si queremos argumentos de entrada mas generales. Por ahora, es una construccion
% particular. 
% - Completar la descripcion de esta funcion. 
%

if nargin == 5
    graficar_malla = 1;
end

if nargin == 6
    truncated = 0;
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
method_data.nsub_coarse= initial_num_el*ones(1,dim);       % Number of subdivisions
method_data.nsub_refine= 2*ones(1,dim);  
method_data.nquad      = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.space_type = space_type;           % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated = truncated;

[hmsh, hspace] = adaptivity_initialize_laplace (problem_data, method_data);

nref = numel(cells);

adaptivity_data.flag = 'elements';

for ref = 1:nref
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MARK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    marked = cell(ref,1);
    marked{ref} = cells{ref};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REFINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    
end

if graficar_malla
    hmsh_plot_cells (hmsh);
end
