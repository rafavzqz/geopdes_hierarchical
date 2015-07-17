function C = split_basis(proj)
%
% function C = split_basis(proj)
%
% Compute the matrix for the change of basis from level lev to level lev+1
%
% Input: 
%           proj = hspace.Proj(lev,:)
%
% Ouput:
%           C: (N_{lev+1} x N_lev matrix), where N_lev is the number of
%           basis functions at level lev. Each column contains the
%           coefficients for writing the corresponding function of level
%           lev as a linear combination of the basis functions of level
%           lev+1
%

dim = numel(proj);

switch dim
    case 1, C = proj{1};
    case 2, C = kron(proj{2}, proj{1});
    case 3, C = kron(proj{3},kron(proj{2},proj{1}));
end
% Provisorio:
%ppp = find(abs(C)<1e-10)
%C(ppp) = 0; 