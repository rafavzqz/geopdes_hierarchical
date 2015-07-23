function [hmsh, hspace] = tp2hier (msh, space, geometry, space_type, boundary)
%
% function [hmsh, hspace] = tp2hier (msh, space, geometry, space_type, boundary)
%
% This function initializes the structures hmsh and hpsace from the
% tensor-product mesh and space. For the initialization of the mesh, the
% function tp2hier_msh is used.
%
% INPUT
%           msh: the coarsest mesh (level 1)
%           space: the spline space for the coarsest space (level 1)
%           geometry: 
%           space_type: 0 (simplified basis), 1 (full basis)
%           boundary: true or false, default: true. (Fill the information
%           for the boundaries).
%
% OUTPUT
%                       struct hmsh (Hierarchical mesh)
% 
%               FIELD_NAME    TYPE               DESCRIPTION
%               ndim           (scalar)          number of parametric directions
%               rdim           (scalar)
%               nlevels       (scalar)           the number of levels
%               nsub          (1 x ndim array)           number of subdivisions on each level
%               mesh_of_level (1 x nlevels mesh) Cartesian mesh of each level (see msh_cartesian)
%               nel           (scalar)           total number of active cells  
%               nel_per_level (1 x nlevels array) number of active cells on each level
%               globnum_active (nel x (dim+1))  global tensor-product numbering of active cells and their corresponding level
%               active        (1 x nlevels cell-array) List of active elements on each level
%               removed       (1 x nlevels cell-array) List of removed cells on each level
%               msh_lev     (nlevels x 1 cell-array) msh_lev{ilev} is a structure
%               geometry
%               boundary    
%
%
%                   hspace: structure for the hierarchical space. It contains the following.
%   
% FIELD_NAME     SIZE                         DESCRIPTION
% ndim           (scalar)           number of parametric directions
% degree        (1 x ndim array)    polynomial degree in each parametric direction
% ndof          (scalar)           total number of active functions 
% nlevels       (scalar)           the number of levels
% space_of_level (1 x nlevels)                tensor product space of each level, with 1d evaluations on the mesh of the same level (see sp_bspline)
% Proj           (hmsh.nlevels-1 x ndim cell-array) the coefficients relating 1D splines of two consecutive levels
%     Proj{l,i} is a matrix of size N_{l+1} x N_l where N_l is the number of univariate functions of level l in the direction l, such that
% A function B_{k,l} = \sum_j c^k_j B_{j,l+1}, and c_j = Proj{l,i}(j,k)
% globnum_active (ndof x (ndim+1))  global tensor-product numbering of active functions and their corresponding level
% ndof_per_level (1 x nlevels array) number of active functions on each level
% active        (1 x nlevels cell-array) List of active functions on each level
% coeff         (ndof x 1)         coefficientes to form the partition of the unity in the hierarchical space
% removed       (1 x nlevels cell-array) List of removed functions on each level
% C             (1 x hmsh.nlevels cell-array) Matrices for changing basis (see compute_matrices_for_changing_basis.m)
% sp_lev        (hmsh.nlevels x 1 cell-array) sp_lev{ilev} is a structure
%

if nargin == 4
    boundary = true;
end

hmsh = tp2hier_msh (msh, geometry, boundary);

hspace = tp2hier_space (hmsh, space, space_type, boundary);

