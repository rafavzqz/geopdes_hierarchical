function hmsh = tp2hier_msh (msh, geometry, boundary)
%
% function hmsh = tp2hier_msh (msh, geometry, boundary)
%
% This function initializes a struct hmsh from the tensor product mesh msh
%
% INPUT
%                       msh: the coarsest mesh (level 1)
%                       geometry: 
%                       boundary: true or false, default: true. (Fill the information for the boundaries of the mesh).
% OUPUT
%                       struct hmsh (Hierarchical mesh)
% 
%               FIELD_NAME    TYPE               DESCRIPTION
%               ndim           (scalar)          number of parametric directions
%               rdim           (scalar)
%               nlevels       (scalar)           the number of levels
%               nsub          (1 x ndim array)       
%               mesh_of_level (1 x nlevels mesh) Cartesian mesh of each level (see msh_cartesian)
%               nel           (scalar)           total number of active cells  
%               nel_per_level (1 x nlevels array) number of active cells on each level
%               globnum_active (nel x (dim+1))  global tensor-product numbering of active cells and their corresponding level
%               active        (1 x nlevels cell-array) List of active elements on each level
%               deactivated       (1 x nlevels cell-array) List of removed cells on each level
%               msh_lev     (nlevels x 1 cell-array) msh_lev{ilev} is a structure
%               geometry
%               boundary    
%
% ATENCION: AQUI nsub == 2 (Dyadic refinement)  <-- modificar
%


% I store nsub, to understand the relation between levels
% nsub = 2 is for dyadic refinement

if nargin == 2
    boundary = true;
end

hmsh.ndim = msh.ndim;
hmsh.rdim = msh.rdim;
hmsh.nlevels = 1;
hmsh.nsub = 2 * ones(1,hmsh.ndim); % Provisorio
hmsh.mesh_of_level = msh;
hmsh.nel = msh.nel;
hmsh.nel_per_level = [msh.nel];

aux = cell(hmsh.ndim,1);
[aux{:}] = ind2sub(msh.nel_dir,1:msh.nel);
hmsh.globnum_active = [ones(msh.nel,1) cell2mat(aux)'];
hmsh.active{1} = (1:msh.nel)';
hmsh.deactivated{1} = zeros(0,1);
hmsh.msh_lev{1} = msh_evaluate_element_list(hmsh.mesh_of_level(1), hmsh.active{1});
hmsh.geometry = geometry;

if (boundary && msh.ndim > 1)
    for iside = 1:2*msh.ndim
        hmsh.boundary(iside) = tp2hier_msh (msh.boundary(iside), geometry.boundary(iside), false);
    end
else
    hmsh.boundary = [];
end
