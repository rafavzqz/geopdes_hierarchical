function [ genealogy ] = get_ancestors(hmsh,Q_ind,lev_Q,lev)
%GET_ANCESTORS(hmsh,Q_ind,lev_Q,lev) returns the indices of the ancestor(s) of level lev of the
%cell(s) on the level lev_Q with indices Q_ind (hierarchical mesh with p-adic refinement)

if lev>=lev_Q
    error('lev>= lev_Q: use get_children')
end
  
if lev==lev_Q-1
    [genealogy,  ~] = hmsh_get_parent (hmsh, lev_Q, Q_ind);
else
    genealogy=get_ancestors(hmsh,hmsh_get_parent(hmsh, lev_Q, Q_ind),lev_Q-1,lev);
end
end