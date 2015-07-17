function value = check_partition_of_the_unity(hmsh, hspace)
%
% function value = check_partition_of_the_unity
%

Z = hspline_eval(hspace.coeff, hmsh, hspace, 0);

value = (max(abs(Z(:)-1)) < 1e-7);