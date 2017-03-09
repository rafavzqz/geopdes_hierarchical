function [ B ] = bzrproj_el_href( p, C_source, C_target, index_target, n_sub )
%BZRPROJ_EL_HREF:  Bézier projection operator (h-refinement)
% This projection operator returns the global h/2 refinement of an
% univariate B-spline element.
%
% Calling Sequence:
%
%   B = bzrproj_el_href( p, C_source, C_target, index_target_elements )
%
%    INPUT:
%      p                        - polynomial degree
%      C_source                 - Bézier extraction operator of the source mesh
%      C_target                 - Bézier extraction operator of the target mesh
%      element.
%      index_target             - starting target index
%      n_sub                    - number of subdivision for refinement
%
%    OUTPUT:
%
%      B - Bézier projection operator [p+1 x p+1 x numel(index_target_elements)]
%
%   Adapted from "Bézier projection: A unified approach for local
%   projection and quadrature-free refinement and coarsening of NURBS and
%   T-splines with particular application to isogeometric design and
%   analysis", D.C. Thomas et al., Comput. Methods Appl. Mech. Engrg. 284
%   (2015)55-105
%
%   Algorithm 4.13 pp. 91-93
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
xi = linspace(-1, 1, n_sub+1);
B = zeros(p+1, p+1, n_sub);

% loop over target elements
for i=1:n_sub
    A = brnsttransf(xi(i), xi(i+1), p);
    
    B(:,:,i) = inv(C_target(:,:,index_target+(i-1)))'*A*C_source(:,:)';
end

end