function [ B ] = brnstpoly_coeff( u, p, i )
%BRNSTPOLY_COEFF:  Bersnstein polynomial coefficients
%
% Calling Sequence:
%
%   B = brnstpoly_coeff( u, p, i )
%
%    INPUT:
%      u  - parametric coordiante in the interval [-1, 1]
%      p  - polynomial degree
%      i  - index of the polynomial
%
%    OUTPUT:
%
%      B - ith Berstein polynomial
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

bin_coeff = factorial(p)/(factorial(i-1)*factorial(p+1-i));
one_u = ones(size(u,1),size(u,2));
B = 1/2^p * bin_coeff * (one_u-u).^(p-(i-1)) * ((one_u+u).^(i-1))';


end


