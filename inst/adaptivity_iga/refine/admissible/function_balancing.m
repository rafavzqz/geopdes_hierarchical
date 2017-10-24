%FUNCTION_BALANCING evaluate a list of elements which have to be further
%refined to obtained a balanced function mesh with function support
%balancing parameter n_fsb.
%
% the function implements Algorithm 7.2 of:
% Lorenzo et al., "Hierarchically refined and coarsened splines for moving
% interface problems, with particular application to phase field models of
% prostate tumor growth", 2017, CMAME pp. 1-46
%
%   balance_element_list = function_balancing( hmsh, hspace, n_fsb )
%
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
function balance_element_list  = function_balancing( hmsh, hspace, n_fsb )

balance_element_list = cell(hmsh.nlevels,1);
% loop over levels
for iLevel = hmsh.nlevels:-1:2
    active =  hmsh.active{iLevel};
    anchestor_level = iLevel - n_fsb - 1;
        % loop over all active element of the ith level
        for el = 1:hmsh.nel_per_level(iLevel)
            if anchestor_level > 0
                active_supp = intersect( hmsh.active{anchestor_level}, support_ext(hmsh, hspace, active(el), iLevel, anchestor_level ) );
                balance_element_list{anchestor_level} = [balance_element_list{anchestor_level}; active_supp]; 
            end
            if ( n_fsb > 1 && iLevel > 2 )
                ring_neighbours = ring_neighbors( hmsh, active(el), iLevel, 1 );
                ring_neighbours = intersect( hmsh.active{iLevel-2}, ring_neighbours);
                balance_element_list{iLevel - 2} = [balance_element_list{iLevel - 2}; ring_neighbours];
            end
        end    
end

end

