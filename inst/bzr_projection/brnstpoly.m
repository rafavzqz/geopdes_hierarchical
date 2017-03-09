function [ B ] = brnstpoly( u, p )
%BRNSTPOLY:  Bersnstein polynomials vector
%
% Calling Sequence:
%
%   B = brnstpoly( u, p )
%
%    INPUT:
%      u  - parametric coordiante in the interval [-1, 1]
%      p  - polynomial degree
%
%    OUTPUT:
%
%      B - Berstein polynomials vector
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

% initialize the vector of p+1 polynomials
B = zeros(1, (p+1));

for i=1:p+1
   bin_coeff = factorial(p)/(factorial(i-1)*factorial(p+1-i));
   b_i_p = 1/2^p * bin_coeff * (1-u)^(p-(i-1)) * (1+u)^(i-1);
   
   % push_back the ith polynomial
   B(i) = b_i_p;     
   
end % end for loop i

end


%! test:
%! isequal(brnstpoly( 0, 2 ), [1/4, 1/2, 1/4] )
%! isequal(brnstpoly( 0.5, 4 ), [1/256, 3/64, 27/128, 27/64, 81/256])
%! isequal(brnstpoly( -0.25, 3 ), [125/512, 225/512, 135/512, 27/512])
%! isequal(brnstpoly( 1, 2 ), [0, 0, 1])
%! isequal(brnstpoly( -1, 5 ), [1, 0, 0, 0, 0, 0])
