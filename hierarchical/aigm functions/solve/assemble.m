function [rhs, K, M] = assemble(hmsh, hspace, f)
%
% function [rhs, K, M] = assemble(hmsh, hspace, f)
%
%

  K = spalloc (hspace.ndof, hspace.ndof, prod(hspace.degree+1)*hspace.ndof);
  M = K;
  rhs = zeros (hspace.ndof, 1);
  
  C = hspace.C;
  msh_lev = hmsh.msh_lev;
  sp_lev = hspace.sp_lev;
  
  ndof_per_level = hspace.ndof_per_level;
  dif = hmsh.nlevels - hspace.nlevels;
  if dif
      ndof_per_level = [ndof_per_level(:); zeros(dif,1)];
  end
  
  ndofs = 0;
  for ilev = 1:hmsh.nlevels % Active levels
      ndofs = ndofs + ndof_per_level(ilev);
      if msh_lev{ilev}.nel
          x = cell(msh_lev{ilev}.rdim,1);
          for idim = 1:msh_lev{ilev}.rdim
              x{idim} = msh_lev{ilev}.geo_map(idim,:,:);
          end
          K_lev = op_gradu_gradv (sp_lev{ilev}, sp_lev{ilev}, msh_lev{ilev}, ones(msh_lev{ilev}.nqn, msh_lev{ilev}.nel));
          M_lev = op_u_v (sp_lev{ilev}, sp_lev{ilev}, msh_lev{ilev}, ones(msh_lev{ilev}.nqn, msh_lev{ilev}.nel));
          b_lev = op_f_v (sp_lev{ilev}, msh_lev{ilev}, f(x{:}));
          
          dofs = 1:ndofs; % Active dofs from level 1 to ilev, in the numbering of the hierarchical space, whatever it is
          % C can also be restricted to the non-vanishing dofs in M_lev
          K(dofs,dofs) = K(dofs,dofs) + C{ilev}'*K_lev*C{ilev};
          M(dofs,dofs) = M(dofs,dofs) + C{ilev}'*M_lev*C{ilev};
          
          rhs(dofs) = rhs(dofs) + C{ilev}'*b_lev;
      end
  end
  