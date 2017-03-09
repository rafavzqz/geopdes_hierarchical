function [ B ] = bzrproj_href(  p, C_source, C_target, n_sub )
%BZRPROJ_HREF:  Bézier projection operator (h-refinement)
% This projection operator returns the global h/2 refinement of an
% univariate B-spline patch.
%
% Calling Sequence:
%
%   B = bzrproj_el_href( p, C_source, C_target )
%
%    INPUT:
%      p                     - polynomial degree
%      C_source              - Bézier extraction operator of the source (coarse) mesh
%      C_target              - Bézier extraction operator of the target (fine) mesh
%      n_sub                 - number of children of each element
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

% evaluate source and target parameters
ndofs_source = size(C_source, 3)+(p-1);         % number of functions of the source mesh
ndofs_target = size(C_target, 3)+(p-1);         % number of functions of the target mesh
nel_source = ndofs_source - p;                  % number of elements of the source mesh

% initialize global Bézier projection operator
B = sparse(ndofs_target, ndofs_source);

% loop over target elements and assembly the corresponding Bézier
% projection operator
for i=1:nel_source
    % index target elements
    start_index = i + (i-1) * (n_sub-1) ;
    % assembly the operators of each refined element
    B_target_el = zeros(n_sub + p, p+1);
    % projection operator of the ith target element
    B_target = bzrproj_el_href(  p, C_source(:,:,i), C_target, start_index, n_sub);
    for j=1:n_sub
        B_target_el(j:j+p,:) = B_target_el(j:j+p,:) + B_target(:,:,j);
    end
    
        % assembly global projection operator 
    B(i + (n_sub-1)*(i-1) : (i-1) + p + n_sub + (n_sub-1)*(i-1), i:i+p) =...
        B(i + (n_sub-1)*(i-1) : (i-1) + p + n_sub + (n_sub-1)*(i-1), i:i+p) + B_target_el;
end

%% Smoothing ==============================================================
%
%     Weights for the Smoothing are tooken from "Convergence of an
%     efficient local least-squares fitting metho for bases with compact
%     support"; Govindjee, S., Strain, J., Mitchell T. J. and Taylor R. L.;
%     Comput. Methods. Appl. Mech. Engrg. 213-216 (2012) 84-92.
%
%

% evaluate weights for smoothing
weights = zeros(ndofs_target, 1);
for i=1:ndofs_target
    if (i < p + 1)
        weights(i) = 1/i;
    elseif (i > ndofs_target - p)
        weights(i) = 1/(ndofs_target-i+1);
    else
        weights(i) = 1/(p+1);
    end
    B(i,:) = weights(i)*B(i,:);
end

end