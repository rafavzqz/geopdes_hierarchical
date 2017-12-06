function [R, hierc_funct_indices, C, L] = mlrcnst(hmsh, hspace, active, deactivated, lev, ind)

[C, hierc_funct_indices] = mlbzrextr(hmsh, hspace, active, lev, ind);
L = trns_ref(hmsh, hspace, active, deactivated, lev, ind);
R = (L*C)\L;

end

