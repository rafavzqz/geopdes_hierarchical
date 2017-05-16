function [ w ] = wgheval( C, fun_index, el_index, fnc_sup_el_index_vec )
%WGHEVAL: Bézier projection weights valuation
% Evaluate the weights to smooth the projected control values.
%
% Calling Sequence:
%
%   w = wgheval( C, fun_index, el_index, fnc_sup_el_index_vec )
%
%    INPUT:
%      C                    - Bézier extraction operator
%      fun_index            - index of the function to smooth
%      el_index             - index of elements of the source mesh to project
%      fnc_sup_el_index_vec - vector of element index support of the function
%
%    OUTPUT:
%
%      w - weights of Bézier projection operator
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


c_num  = 0.0;
c_den = 0.0;

% loop over control values...
for i=1:size(C,2)
    c_num = c_num + C(fun_index(fnc_sup_el_index_vec==el_index), i, el_index);
end
% ... loop over each support element
for i_el=1:numel(fnc_sup_el_index_vec)
    c_den_i = 0.0;
    % ... and loop again over the control values
    for i=1:size(C,2)
        c_den_i = c_den_i + C(fun_index(i_el), i, fnc_sup_el_index_vec(i_el));
    end
    c_den = c_den_i + c_den;
end

if (c_den==0)
   disp('ERROR: Weight is not a number !!!') 
elseif (c_den==0 && c_num==0)
   disp('ERROR: Weight is zero !!!') 
end
w = c_num/c_den;

end

%! test: 
%! p = 2;
%! knots = [0 0 0 1/3 2/3 1 1 1];
%! C = bzrextr(knots, p);
%! if (isequal(wgheval( C, 1, 1, 1 ), 1) && ...
%!      isequal(wgheval(  C, [2 1], 1, [1 2] ),[3/4]) && ...
%!       isequal(wgheval(  C, [3 2 1], 2, [1 2 3] ),[2/3]) && ...
%!        isequal(wgheval(  C, [1 2], 2, [2 3] ),[1/4]) && ...
%!         isequal(wgheval(  C, 3, 3, 3 ),1))
%!     disp('Sucess !!!');
%! else
%!     disp('Error !!!');
%! end

