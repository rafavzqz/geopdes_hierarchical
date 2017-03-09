function [ C ] = bzrextr( XI, poly )
%BZREXTR:  Bèzier extraction operator
%
% Calling Sequence:
%
%   C = bzrextr( XI, poly )
%
%    INPUT:
%      XI    - knot sequence
%      poly  - polynomial degree
%
%    OUTPUT:
%
%      C - Bézier extraction operator
%
%    Adapted from "Isogeometric Finite Element Analysis based on Bézier 
%    Extraction of NURBS and T-Splines", T. N. Nguyen, Master Thesis, NTNU.
%
%
% Copyright (C) 2017 Davide D'Angella
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

% Initializations
m=length(XI);
a=poly+1;
b=a+1;
nb=1;
C(:,:,1)=eye(poly+1);

% for k=1:m-a+1
% C(:,:,k)=eye(poly+1);
% end
% return;

while b<m
    C(:,:,nb+1)=eye(poly+1);
    i=b;
    % Initialize the next extraction operator
    % Count multiplicity of the knot at location b
    while b<m && XI(b+1)==XI(b)
        b=b+1;
    end
    mult=b-i+1;
    if mult<poly
        numer=XI(b)-XI(a);
        for j=poly:-1:mult+1
            alphas(j-mult)=numer/(XI(a+j)-XI(a));
        end
        r=poly-mult;
        % Update the matrix coefficient for r new knots
        for j=1:r
            save=r-j+1;
            s=mult+j;
            for k=poly+1:-1:s+1
                alpha=alphas(k-s);
                % Form extraction operator
                C(:,k,nb)=alpha*C(:,k,nb)+(1-alpha)*C(:,k-1,nb);
            end
            if b<m
                % Update overlapping coefficients of the next operator
                C(save:j+save,save,nb+1)=C(poly-j+1:poly+1,poly+1,nb);
            end
        end
        % Finished with the current operator.
        % Update indices for the next operator.
        nb=nb+1;
        if b<m
            a=b;
            b=b+1;
        end
    elseif mult==poly
        % In case multiplicity of knot is already p,
        % update indices for the next operator.
        nb=nb+1;
        if b<m
            a=b;
            b=b+1;
        end
    end
end
end

