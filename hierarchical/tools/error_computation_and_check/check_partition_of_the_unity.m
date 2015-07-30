function value = check_partition_of_the_unity(hmsh, hspace)
%
% function value = check_partition_of_the_unity(hmsh, hspace)
%
% This function checks the property of partition of unity on the quadrature
% nodes
%
%

% Z = hspline_eval_old (hspace.coeff, hmsh, hspace, 0);
Z = hsp_eval_hmsh (hspace.coeff, hspace, hmsh);

value = (max(abs(Z(:)-1)) < 1e-7);