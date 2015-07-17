function [hmsh, hspace] = build_hspace_from_cells(dim, p, initial_num_el, cells,flag_whole_basis, boundary, graficar_malla)

if nargin == 6
    graficar_malla = 1;
end

switch dim
    case 1, geometry = geo_load (nrbline ([0 0], [1 0]));
    case 2, geometry = geo_load ('geo_square.txt'); %geometry=geo_load ('geo_ring.txt'); 
    case 3, geometry = geo_load ('geo_cube.txt');
end

[tp_msh, tp_space] = get_initial_msh_and_space(dim, p, initial_num_el,geometry);
[hmsh, hspace] = tp2hier (tp_msh, tp_space, geometry);

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
    
    [hmsh, hspace] = refine (hmsh, hspace, marked, 'elements', flag_whole_basis, boundary);
    
end

if graficar_malla
plot_msh(hmsh)
end
