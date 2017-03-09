function [G] = grmn( p )
%GRMN:  Gramian of Bersnstein basis
%
% Calling Sequence:
%
%   G = grmn( p )
%
%    INPUT:
%      p        - polynomial degree
%
%    OUTPUT:
%
%      G - Gramian matrix of Berstein basis
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

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% initialize the matrix G
G = zeros((p+1), (p+1));

for k=1:p+1
    for j=1:p+1
        G(j,k) = (2/(2*p+1))*(factorial(2*p)/(factorial(j+k-2)*factorial(2*p-j-k+2)))^(-1)*...
            (factorial(p)/(factorial(j-1)*factorial(p-j+1)))*(factorial(p)/(factorial(k-1)*factorial(p-k+1)));
    end % end loop j
end % end for loop k

end

%! test: 
%! p = 2
%! G_ex = [[2/5 1/5 1/15]; [1/5 4/15 1/5]; [1/15 1/5 2/5]];
%! G = grmn( p );
%! if (G == G_ex)
%!     disp('Sucess !!!');
%! else
%!     dips('Error !!!');
%! end

