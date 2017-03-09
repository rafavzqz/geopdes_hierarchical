function [ B ] = bzrproj_el_hcoarse( p, C_source, R_target, h_source, h_target, index_source_element_vec )
%BZRPROJ_EL_HCOARSE:  Bézier projection operator (h-coarsening)
% This projection operator project n elements of a source mesh onto
% an element of a target mesh. The resulting control values have to be 
% smoothed by means of the weights.
%
% Calling Sequence:
%
%   B = bzrproj_hcoarse( p, C_source, R_target, h_source, h_target, taget_el_bound_vec, index_source_element_vec )
%
%    INPUT:
%      p                        - polynomial degree
%      C_source                 - Bézier extraction operator of the source mesh
%      R_target                 - spline reconstruction operator of the target element (inverse of extraction operator)
%      h_source                 - element size of the source mesh
%      h_target                 - element size of the target mesh
%      index_source_element_vec - vector of index source elements
%
%    OUTPUT:
%
%      B - Bézier projection operator [p+1 x p+1 x n+1]
%
%   Adapted from "Bézier projection: A unified approach for local
%   projection and quadrature-free refinement and coarsening of NURBS and
%   T-splines with particular application to isogeometric design and
%   analysis", D.C. Thomas et al., Comput. Methods Appl. Mech. Engrg. 284
%   (2015)55-105
%
%   Algorithm 4.14 pp. 94-95
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

% initialize Bézier projection operator
n = numel(index_source_element_vec);
B = zeros(p+1, p+1, n);

% build up matrix to convert Bézier coefficients of the source Bernstein
% basis into coefficients of the target Bernstein basis

% Gramian matrix and inverse
G = grmn(p);
G_inv = grmninv(p);
% source elements intervals in parent coordinates
xi = linspace(-1, 1, n+1);
% elements volume rate
phi = h_source/h_target;

% loop over elements of the source mesh
for i=1:n
    % transformation matrix from [-1, 1] to taget element boundaries
    A = brnsttransf(xi(i), xi(i+1), p);
    % projection operator
    B(:,:,i) = R_target(:,:)' * phi*G_inv*A'*G * C_source(:,:,index_source_element_vec(i))';
end

end