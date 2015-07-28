function C = compute_matrices_for_changing_basis(lev, ind, ndof, proj)
%
% function C = compute_matrices_for_changing_basis(lev, ind, ndof, proj)
%
% Computes the matrices for the change of basis from the hierarchical basis
% functions up to level lev to the tensor-product basis functions of level
% lev
%
% Input: 
%           lev: level (scalar)
%           ind: (1 x lev cell-array), ind{l} is a vector containing
%           the indices of active functions of level l, for l = 1:lev
%           ndof(l) = hspace.space_of_level(l).ndof, for l = 1:lev
%           proj(l,:) = hspace.Proj(l,:), for l = 1:lev-1
%
% Ouput:
%           C{l}: (N_{lev} x N matrix), where N_l is the number of
%           basis functions at level l and N is the number of active functions up to level l, for l = 1:lev. Each column contains the
%           coefficients for writing the corresponding active function as a linear combination of the basis functions of level
%           l
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% REMARK: We have to put it into hspace class
%


C = cell(lev,1);

% ndof = zeros(lev,1);
% for l = 1:lev
%     ndof(l) = hspace.space_of_level(l).ndof;
% end

C{1} = speye(ndof(1));
C{1} = C{1}(:,ind{1});

for l = 2:lev
    I = speye(ndof(l));
    C{l} = [split_basis(proj(l-1,:))*C{l-1}, I(:,ind{l})];
end