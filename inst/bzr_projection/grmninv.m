function [G_inv] = grmninv( p )
%GRMNINV:  inverse of Gramian of Bersnstein basis
%
% Calling Sequence:
%
%   G_inv = grmninv( p )
%
%    INPUT:
%      p        - polynomial degree
%
%    OUTPUT:
%
%      G_inv - inverse Gramian matrix of Berstein basis
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

% initialize the matrix G_inv
G_inv = zeros((p+1), (p+1));

for k=1:p+1
    for j=1:p+1
        sum_i = 0.0;
        for i=1:min([j k])
            sum_i=(2*i-1)*(factorial(p-i+1)/(factorial(p-j+1)*factorial(p-i+1-p+j-1)))*...
                (factorial(p-i+1)/(factorial(p-k+1)*factorial(p-i+1-p+k-1)))*(factorial(p+i)/(factorial(p-j+1)*factorial(p+i-p+j-1)))*...
                (factorial(p+i)/(factorial(p-k+1)*factorial(p+i-p+k-1))) + sum_i;            
        end
        G_inv(j,k) = (-1)^(j+k)/2 * ((factorial(p)/(factorial(j-1)*factorial(p-j+1)))*(factorial(p)/(factorial(k-1)*factorial(p-k+1))))^(-1)*sum_i;
    end % end loop j
end % end for loop k

end

%! test: 
%! p = 2
%! G_ex = [[2/5 1/5 1/15]; [1/5 4/15 1/5]; [1/15 1/5 2/5]];
%! G_inv = grmninv( p );
%! if (eq(G_inv, inv(G_ex)))
%!     disp('Sucess !!!');
%! else
%!     disp('Error !!!');
%! end

