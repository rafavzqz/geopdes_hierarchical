function [ A ] = brnsttransf( a_tilde, b_tilde, p )
%BRNSTTRANSF:  Bersnstein transformation matrix
% Map the Bernstein polynomials defined over a an interval [a, b] onto a
% new interval [a_tilde, b_tilde].
%
% Calling Sequence:
%
%   A = brnsttransf( a_tilde, b_tilde, p )
%
%    INPUT:
%      a_tilde  - left end of the new domain
%      b_tilde  - right end of the new domain
%      p        - polynomial degree
%
%    OUTPUT:
%
%      A - Berstein transformation matrix
%
%   Adapted from "On the Numerical Condition of Berstein-Bézier Subdivision
%   Processes", R.T. Farouki, C.A. Neff, Mathematics of Computation,
%   Vol.55, No. 192(1990), pp. 637-647.
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

% initialize the matrix A
A = zeros((p+1), (p+1));

for k=1:p+1
    for j=1:p+1
        h = j+k-(p+1);
        A(j,k) = 0.0;         % initial value of the jk entry of A
        for i=max([1 h]):min([j k])
            % evaluate Benstein polynomials vector in a_tilde
            b_a = brnstpoly(a_tilde, p-j+1);
            % evaluate Benstein polynomials vector in b_tilde
            b_b = brnstpoly(b_tilde, j-1);
            
            A(j,k) = A(j,k) + b_b(i)*b_a(k-i+1);
            
        end % end loop i
    end % end loop j
end % end for loop k

end

%! test 1:
%! p = 2
%! a = -1;
%! b =  3;
%! A1_inv_ex = [[1 0 0]; [-1 2 0]; [1 -4 4]];
%! A1_inv = brnsttransf( a, b, p );
%! if (A1_inv == A1_inv_ex)
%!     disp('Sucess !!!');
%! else
%!     disp('Error !!!');
%! end


%! test 2:
%! p = 2
%! a = -1;
%! b =  0;
%! A1_inv_ex = [[1 0 0]; [1/2 1/2 0]; [1/4 1/2 1/4]];
%! A1_inv = brnsttransf( a, b, p );
%! if (A1_inv == A1_inv_ex)
%!     disp('Sucess !!!');
%! else
%!     disp('Error !!!');
%! end


