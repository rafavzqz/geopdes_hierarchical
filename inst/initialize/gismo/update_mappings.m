% UPDATE_MAPPINGS: Given a geometry (obtained with geo_load), updates a
% hierarchical or tensor product mesh with the geometry mapping.
%
% function msh = update_mappings (msh, geometry, rdim)
%
% INPUT: 
%   msh: msh_cartesian or hierarchical_mesh object
%   geometry: geometry structure obtained from geo_load (see geo_load.m)
%   rdim (optional): dimension of the physical space. Default: geometry.rdim
%
% OUTPUT: 
%   msh: msh_cartesian or hierarchical_mesh object with updated mappings.
%
% Copyright (C) 2018 Ondine Chanon


function msh = update_mappings (msh, geometry, rdim)

if nargin < 3
    rdim = geometry.rdim;
end
if isa (msh, 'msh_cartesian')
    msh.rdim = rdim;
    msh.map = geometry.map;
    msh.map_der = geometry.map_der;
    if isfield (geometry, 'map_der2')
        msh.map_der2 = geometry.map_der2;
    end

    if ~isempty (msh.boundary)
        assert (length (geometry.boundary) == length (msh.boundary))

        for ibd = 1:length (geometry.boundary)
            msh.boundary(ibd).rdim = rdim;
            if msh.boundary(ibd).ndim > 0
                msh.boundary(ibd).map = geometry.boundary(ibd).map;
                msh.boundary(ibd).map_der = geometry.boundary(ibd).map_der;

                if isfield (geometry.boundary(ibd), 'map_der2')
                    msh.boundary(ibd).map_der2 = geometry.boundary(ibd).map_der2;
                end
            end
        end
    end
    
elseif isa (msh, 'hierarchical_mesh')
    msh.rdim = rdim;
    for ilev = 1:msh.nlevels
        msh.mesh_of_level(ilev).rdim = rdim;
        msh.mesh_of_level(ilev).map = geometry.map;
        msh.mesh_of_level(ilev).map_der = geometry.map_der;
        msh.mesh_of_level(ilev).map_der2 = geometry.map_der2;
        msh.msh_lev{ilev} = msh_evaluate_element_list ...
            (msh.mesh_of_level(ilev), msh.active{ilev});
    end
    for ibd = 1:length (msh.boundary)
        msh.boundary(ibd).rdim = rdim;
        if msh.boundary(ibd).ndim > 0
            for ilev = 1:msh.boundary(ibd).nlevels
                msh.boundary(ibd).mesh_of_level(ilev).rdim = rdim;
                msh.boundary(ibd).mesh_of_level(ilev).map = ...
                    geometry.boundary(ibd).map;
                msh.boundary(ibd).mesh_of_level(ilev).map_der = ...
                    geometry.boundary(ibd).map_der;
                if isfield(geometry.boundary(ibd), 'map_der2')
                    msh.boundary(ibd).mesh_of_level(ilev).map_der2 = ...
                        geometry.boundary(ibd).map_der2;
                end
                msh.boundary(ibd).msh_lev{ilev} = msh_evaluate_element_list ...
                    (msh.boundary(ibd).mesh_of_level(ilev), msh.boundary(ibd).active{ilev});
            end
        end
    end
        
end

end