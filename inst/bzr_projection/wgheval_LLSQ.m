function weight = wgheval_LLSQ( funs_support )
%WGHEVAL_LLSQ: Bézier projection weights valuation
% Evaluate the weights to smooth the projected control values.
%
% Calling Sequence:
%
%   weight = wgheval_LLSQ( p, ndofs, index_dof )
%
%    INPUT:
%
%      funs_support        - cell array of elements supporting the function
%
%    OUTPUT:
%
%      weight - weights of Bézier projection operator
%
%   Adapted from "Convergence of an
%     efficient local least-squares fitting metho for bases with compact
%     support"; Govindjee, S., Strain, J., Mitchell T. J. and Taylor R. L.;
%     Comput. Methods. Appl. Mech. Engrg. 213-216 (2012) 84-92.
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

weight = 1/numel(funs_support) ;

end

%! test: 
%! if (isequal(wgheval_LLSQ( 2, 5, 5 ), 1) && ...
%!      isequal(wgheval_LLSQ( 4, 10, 5 ), 1/5) && ...
%!       isequal(wgheval_LLSQ( 2, 4, 3 ), 1/2) && ...
%!        isequal(wgheval_LLSQ( 3, 6, 4 ),1/3))
%!     disp('Sucess !!!');
%! else
%!     disp('Error !!!');
%! end

