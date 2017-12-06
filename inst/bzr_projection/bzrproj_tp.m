function [ B ] = bzrproj_tp( p, hmsh, source_lev, target_lev, child_subIndices, base_subIndices )
%BZRPROJ_EL_HCOARSE:  Bezier projection operator (h-coarsening)
% This projection operator project n elements of a source mesh onto
% an element of a target mesh. The resulting control values have to be 
% smoothed by means of the weights.
%
% Calling Sequence:
%
%   B = bzrproj_hcoarse( p, hmshc )
%
%    INPUT:
%      p                        - polynomial degree
%      hmsh                     - hierarchical mesh structure
%      lev                      - level
%
%    OUTPUT:
%
%      B - Bezier projection operator 
%
%   Adapted from "Bezier projection: A unified approach for local
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

    
% tensor product 
B = 1;

for iDir=1:hmsh.ndim
    scaling = hmsh.nsub(iDir);
    % elements volume rate
    phi = 1/scaling^(source_lev-target_lev);
    index = child_subIndices(iDir) - (base_subIndices(iDir)-1)*scaling^(source_lev-target_lev);
    % Gramian matrix and inverse
    G = grmn(p(iDir));
    G_inv = grmninv(p(iDir));
    % source elements intervals in parent coordinates
    xi = linspace(-1, 1,hmsh.nsub(iDir)^(source_lev-target_lev)+1 ); % for bisection only
    % transformation matrix from [-1, 1] to taget element boundaries
    A = brnsttransf(xi(index), xi(index+1), p(iDir));
    % projection operator
    B_dir = phi*G_inv*A'*G ;
    B = kron(B_dir,B);
end

end