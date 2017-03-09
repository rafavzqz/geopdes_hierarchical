function B = bzrproj_hcoarse( p, C_source, C_target, n_sub )
%BZRPROJ_HCOARSE: global Bézier projection operator (h-coarsening)
% This projection operator projects the control values of the source (fine) 
% mesh onto the ones of the target (coarse) mesh.
%
% Calling Sequence:
%
%   B = bzrproj_hcoarse( p, C_source, C_target, taget_el_bound_vec )
%
%    INPUT:
%      p                        - polynomial degree
%      C_source                 - Bézier extraction operator of the source mesh
%      C_target                 - Bézier extraction operator of the source mesh
%      n_sub                    - number of children per element
%
%    OUTPUT:
%
%      B - global Bézier projection operator [ndof_target x ndof_source]
%
%   The function is implemented for the univariate case and uniform knot spans.
%   Smoothing coeffients are taken as: w_A = 1/n_A_el, where n_A_el is the number
%   of supporting elements of the Ath function.
%
%   Adapted from "Bézier projection: A unified approach for local
%   projection and quadrature-free refinement and coarsening of NURBS and
%   T-splines with particular application to isogeometric design and
%   analysis", D.C. Thomas et al., Comput. Methods Appl. Mech. Engrg. 284
%   (2015)55-105
%
% Copyright (C) 2017 Massimo Carraturo
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
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
nel_target = ndofs_target - p;                  % number of elements of the target mesh
h_source = 1/nel_source;                        % element size of the source mesh
h_target = 1/nel_target;                        % element size of the target mesh


% initialize global Bézier projection operator
B = sparse(ndofs_target, ndofs_source);  

% index source elements vector
source_dofs_to_coarse = (ndofs_source - ndofs_target) * n_sub + p;
index_source_element_vec = linspace(1,source_dofs_to_coarse,source_dofs_to_coarse);

% build up matrix to convert Bézier coefficients of the source Bernstein
% basis into coefficients of the target Bernstein basis

% loop over target elements and assembly the corresponding Bézier
% projection operator
for i=1:(nel_source-nel_target)
    % index source elements
    start_index = i + (i-1) * (n_sub-1) ;
    end_index = start_index + n_sub - 1;
    % projection operator of the ith target element
    B_target = bzrproj_el_hcoarse( p, C_source, inv(C_target(:,:,i)), h_source, h_target,...
        index_source_element_vec(start_index:end_index));
    % assembly the operators of each sub-element 
    B_target_el = zeros(p+1, n_sub + p);
    for j=1:n_sub
        B_target_el(:,j:j+p) = B_target_el(:,j:j+p) + B_target(:,:,j);
    end
    % assembly global projection operator 
    B(i:i+p, i + (n_sub-1)*(i-1) : (i-1) + p + n_sub + (n_sub-1)*(i-1)) =...
        B(i:i+p, i + (n_sub-1)*(i-1) : (i-1) + p + n_sub + (n_sub-1)*(i-1)) + B_target_el;
end

%% Smoothing ==============================================================
%
%     Weights for the Smoothing are tooken from "Convergence of an
%     efficient local least-squares fitting method for bases with compact
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
