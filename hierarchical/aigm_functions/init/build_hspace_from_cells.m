function [hmsh, hspace] = build_hspace_from_cells(dim, p, initial_num_el, cells,space_type, graficar_malla)
%
% function [hmsh, hspace] = build_hspace_from_cells(dim, p, initial_num_el, cells, space_type, graficar_malla)
%
% This function fills hmsh and hspace. The active cells are given in
% cells{lev}, for lev = 1,2,...
%
% ATENCION: Completar la descripcion de esta funcion y decidir si queremos argumentos de entrada mas generales. Por ahora, es una construccion
% particular
%

if nargin == 5
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
method_data.nsub_coarse= initial_num_el*ones(1,dim);       % Number of subdivisions
method_data.nsub_refine= 2*ones(1,dim);  
method_data.nquad      = method_data.degree+1;       % Points for the Gaussian quadrature rule
method_data.space_type = space_type;           % 0: , 1: Full basis (B-splines)

[hmsh, hspace] = adaptivity_initialize (problem_data, method_data);

nref = numel(cells);

adaptivity_data.flag = 'elements';

for ref = 1:nref
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MARK
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    marked = cell(ref,1);
    aux = cell (hmsh.ndim, 1);
    for idim = 1:hmsh.ndim
        if ~isempty(cells{ref})
            aux{idim} = cells{ref}(:,idim);
        else
            aux{idim} = [];
        end
    end
    marked{ref} = sub2ind (hmsh.mesh_of_level(ref).nel_dir, aux{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% REFINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [hmsh, hspace] = adaptivity_refine (hmsh, hspace, marked, adaptivity_data);
    
end

if graficar_malla
    hmsh_plot_cells (hmsh, 1);
end
