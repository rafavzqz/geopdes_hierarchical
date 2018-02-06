% MATRIX_BASIS_CHANGE__: compute the subdivision matrix between two consecutive levels.
%        This method is intended to remain private.
%
% function C = matrix_basis_change__ (hspace, lev)
%
% Compute the new matrices to represent functions of level "lev-1"
% as linear combinations of splines (active and inactive) of level "lev"
%
% INPUT:  
%
%   hspace: an object of the class hierarchical_space_mp
%   lev:    the level for which we compute the matrix
%
% OUTPUT:
%
%   C:    matrix to change basis from level lev-1 to level lev
%
% Copyright (C) 2015, 2016 Eduardo M. Garau, Rafael Vazquez
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

function C = matrix_basis_change__ (hspace, lev, ind_coarse)

% We need to change this function, to compute the coefficients in Mario's notes (refinement masks)
% It will write C^1 basis function of level lev-1, as linear combinations of C^1 basis functions of level lev.
% The output matrix must have size: ndof_{l} x ndof_{l-1}
% The coefficients are computed using the projection matrices for
%  univariate spaces in hspace.Proj, hspace.Proj0, hspace.Proj1
%    Proj  (hmsh.nlevels-1 x npatch cell-array) 
%           the coefficients relating 1D splines of two consecutive levels for each patch
%           Proj{l,i} is a cell-array of dimension ndim, with the information for
%           the univariate Projectors on the patch (see also hierarchical_space)
%    Proj0  Like Proj, for univariate splines of degree p, regularity r+1
%    Proj1  Like Proj, for univariate splines of degree p-1, regularity r
% For the numbering at each level, we can't use gnum anymore, as it does not exist for the C1 spaces
% Instead, we have to retrieve this information from hspace.Cpatch, which
%  writes the C^1 basis functions as linear combinations of standard B-splines in the patch


if (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_scalar'))
  is_scalar = true;
  ndim = size (hspace.Proj{1}, 2);
% elseif (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_vector'))
%   is_scalar = false;
%   ndim = size (hspace.Proj{1}, 2);
else
  error ('Unknown space type')
end

npatch = hspace.space_of_level(1).npatch;
ndim = size (hspace.Proj{1}, 2);

% For now we only use B-splines, and scalar spaces

% We have to change all this part following Mario's notes, and using Proj, Proj0 and Proj1
% How are the three types of basis functions organized (indexed) in the hspace object?
if (nargin < 3)
    
  C = sparse (hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof); %zero sparse matrix 

  %Cesare's note: this defines the dimension of the matrix C (considers full multipatch spaces on consecutive levels)
  % in the C^0 case, for each patch you compute the classical matrix going
  % from level l-1 to level l (tensor product basis), and then put it in
  % the right palce in the big matrix, according to the global numbering of
  % the dofs of the patch.
  % What does Proj (and Proj0, Proj1) contain in our case?
  for ipatch = 1:npatch  
      
    %first we construct the matrices containing the coefficients lambda, mu, nu (usual refinement coefficients for B-splines)  
    Lambda{ipatch} = 1;
    Proj = hspace.Proj{lev-1, ipatch};
    for idim = 1:ndim
      Lambda{ipatch} = kron (Proj{idim}, Lambda{ipatch});
    end
%     Mu{ipatch} = 1;
    Proj0 = hspace.Proj0{lev-1, ipatch};
%     for idim = 1:ndim
%       Mu{ipatch} = kron (Proj0{idim}, Mu{ipatch});
%     end   
%     Nu{ipatch} = 1;
    Proj1 = hspace.Proj1{lev-1, ipatch};
%     for idim = 1:ndim
%       Nu{ipatch} = kron (Proj1{idim}, Nu{ipatch});
%     end    
    
    %get interf_dir=1,2 according to the interface being horizontal or vertical
    all_dir=1:ndim; % possible directions
    side = hspace.space_of_level(lev-1).interfaces.side1;
    interf_dir_orthogonal = ceil (side/2);  %direction orthogonal to the interface
    interf_dir_parallel = setdiff (all_dir, interf_dir_orthogonal);  %direction(s) parallel to the interface
    interf_dir = interf_dir_parallel;
    degree=hspace.space_of_level(lev-1).sp_patch{ipatch}.degree(interf_dir); %degree
    
    %define the 1D B-spline space parallel to the interface: spn (p,r)
    %(actually it is enough to know the number of degrees of freedom of this space)
    ndof_dir_spn=hspace.space_of_level(lev-1).sp_patch{ipatch}.ndof_dir;
    ndof_spn=ndof_dir_spn(interf_dir);  
    ndof_Bsp=prod(ndof_dir_spn);  %dimension of bivariate space
    ndof_dir_spn_ref=hspace.space_of_level(lev).sp_patch{ipatch}.ndof_dir; %same for the finer level 
    ndof_spn_ref=ndof_dir_spn_ref(interf_dir);     %same for the finer level
    ndof_Bsp_ref=prod(ndof_dir_spn_ref);   %same for the finer level
    
    %get the indices ind0 and ind1, corresponding, in the matrix CC (see sp_multipatch_C1), 
    %to the B-splines of patch iptc used in the
    %sums of (12) and (13) of Mario's notes (use sub2ind with spn)
    %(needs to be fixed and checked)
      if (side == 1)
        ind0 = sub2ind (ndof_dir_spn, ones(1,ndof_spn), 1:ndof_spn);
        ind1 = sub2ind (ndof_dir_spn, 2*ones(1,ndof_spn), 1:ndof_spn);
        ind0_ref = sub2ind (ndof_dir_spn_ref, ones(1,ndof_spn_ref), 1:ndof_spn_ref);  %same for the finer level
        ind1_ref = sub2ind (ndof_dir_spn_ref, 2*ones(1,ndof_spn_ref), 1:ndof_spn_ref); %same for the finer level       
      elseif (side == 2)
        ind0 = sub2ind (ndof_dir_spn, ndof_dir_spn(1) * ones(1,ndof_spn), 1:ndof_spn);
        ind1 = sub2ind (ndof_dir_spn, (ndof_dir_spn(1)-1) * ones(1,ndof_spn), 1:ndof_spn);
        ind0_ref = sub2ind (ndof_dir_spn_ref, ndof_dir_spn_ref(1) * ones(1,ndof_spn_ref), 1:ndof_spn_ref);  %same for the finer level
        ind1_ref = sub2ind (ndof_dir_spn_ref, (ndof_dir_spn_ref(1)-1) * ones(1,ndof_spn_ref), 1:ndof_spn_ref); %same for the finer level       
      elseif (side == 3)
        ind0 = sub2ind (ndof_dir_spn, 1:ndof_spn, ones(1,ndof_spn));
        ind1 = sub2ind (ndof_dir_spn, 1:ndof_spn, 2*ones(1,ndof_spn));
        ind0_ref = sub2ind (ndof_dir_spn_ref, 1:ndof_spn_ref, ones(1,ndof_spn_ref));  %same for the finer level
        ind1_ref = sub2ind (ndof_dir_spn_ref, 1:ndof_spn_ref, 2*ones(1,ndof_spn_ref));  %same for the finer level        
      elseif (side == 4)
        ind0 = sub2ind (ndof_dir_spn, 1:ndof_spn, ndof_dir_spn(2) * ones(1,ndof_spn));
        ind1 = sub2ind (ndof_dir_spn, 1:ndof_spn, (ndof_dir_spn(2)-1) * ones(1,ndof_spn));
        ind0_ref = sub2ind (ndof_dir_spn_ref, 1:ndof_spn_ref, ndof_dir_spn_ref(2) * ones(1,ndof_spn_ref));  %same for the finer level
        ind1_ref = sub2ind (ndof_dir_spn_ref, 1:ndof_spn_ref, (ndof_dir_spn_ref(2)-1) * ones(1,ndof_spn_ref));  %same for the finer level       
      end      
    
    %define the 1D B-spline spaces parallel to the interface: sp0 (p,r+1) and sp1 (p-1,r)
    %(actually it is enough to know the number of degrees of freedom of these two spaces)
    ndof_0_C1=length(hspace.space_of_level(lev-1).knots0_patches{ipatch}{interf_dir})-degree-1;
    %ndof_1_C1=length(hspace.space_of_level(lev-1).knots1_patches{ipatch}{interf_dir})-degree-1;
    ndof_0_C1_ref=length(hspace.space_of_level(lev).knots0_patches{ipatch}{interf_dir})-degree-1; %same for the finer level
    %ndof_1_C1_ref=length(hspace.space_of_level(lev).knots1_patches{ipatch}{interf_dir})-degree-1; %same for the finer level   
    
    %get the number of internal knots points (of spn):
    %k=ndof_spn-degree-1;   %must be divide by (degree-r)
    
    %REMARK: in CC{iptc}, columns 1:sp0.ndof correspond to \phi_{0,i} and
    %sp0.ndof+1:sp0.ndof+sp1.ndof correspond to \phi_{1,i} 
    %(functions on the interface, which live and have same coefficients in both patches)
    
    %in Cpatch{iptc}, columns 1:sp.ndof_interior represent the regular
    %(interior) B-splines (all patches),
    %sp.ndof_interior+1:sp.ndof_interior+sp0.ndof represent \phi_{0,i}
    %sp.ndof_interior+sp0.ndof+1:sp.ndof represent \phi_{1,i}
    ndof_interior_C1=hspace.space_of_level(lev-1).ndof_interior;
    ndof_interior_C1_ref=hspace.space_of_level(lev).ndof_interior; %same for the finer level
    
    Cpatch=hspace.space_of_level(lev-1).Cpatch;
    
    %in Cpatch{iptc}, rows interior_dofs_per_patch{itpc} represent the interior B-spline
    %ndof_interior=hspace.space_of_level(lev-1).interior_dofs_per_patch{ipatch};
    %ndof_interior_ref=hspace.space_of_level(lev).interior_dofs_per_patch{ipatch}; %same for the finer level
    %ind0 correspond to N_0^(p,r) in (12) and ind1 to N_1^(p,r) in (12)-(13) 
    
%     A0=Cpatch{ipatch}(ind0,(ndof_interior_C1+1):(ndof_interior_C1+ndof_0_C1));% contains the coefficients a_{i,0,j}
%     A1=Cpatch{ipatch}(ind1,(ndof_interior_C1+1):(ndof_interior_C1+ndof_0_C1));% contains the coefficients a_{i,1,j}
%     hA1=Cpatch{ipatch}(ind1,(ndof_interior_C1+ndof_0_C1+1):(hspace.space_of_level(lev-1).ndof));% contains the coefficients \hat a_{i,1,j}
%     
%     zeta{ipatch}=A0;
%     eta{ipatch}=degree*(k+1)*(A0-A1);
%     theta{ipatch}=degree*(k+1)*hA1;

   Aux=Lambda{ipatch}*Cpatch{ipatch}; %this matrix expresses C1 basis functions of lev-1 in terms of B-splines of level lev
    
    %Taking care of the parts corresponding to Lemma 1 (Mario's notes on refinement) 
    ind_start_int_dofs_patch=0; %computing the (starting) global index for the interior dofs of the patch of level lev-1
    ind_start_int_dofs_patch_ref=0; %same for the finer level (lev-1)
    for i=1:ipatch-1;
        ind_start_int_dofs_patch=ind_start_int_dofs_patch+numel(hspace.space_of_level(lev-1).interior_dofs_per_patch{i});
        ind_start_int_dofs_patch_ref=ind_start_int_dofs_patch_ref+numel(hspace.space_of_level(lev).interior_dofs_per_patch{i});
    end
    ind_start_int_dofs_patch=ind_start_int_dofs_patch+1;
    ind_start_int_dofs_patch_ref=ind_start_int_dofs_patch_ref+1;
    ind_end_int_dofs_patch=ind_start_int_dofs_patch+numel(hspace.space_of_level(lev-1).interior_dofs_per_patch{ipatch})-1;
    ind_end_int_dofs_patch_ref=ind_start_int_dofs_patch_ref+numel(hspace.space_of_level(lev).interior_dofs_per_patch{ipatch})-1;

    C(ind_start_int_dofs_patch_ref:ind_end_int_dofs_patch_ref,ind_start_int_dofs_patch:ind_end_int_dofs_patch)=...
        Lambda{ipatch}(setdiff(1:ndof_Bsp_ref,union(ind0_ref,ind1_ref)),setdiff(1:ndof_Bsp,union(ind0,ind1)));  %not all Lambda
    
   %Taking care of the parts corresponding to Lemma 2 (Mario's notes on refinement)
   C(ndof_interior_C1_ref+1:ndof_interior_C1_ref+ndof_0_C1_ref,ndof_interior_C1+1:ndof_interior_C1+ndof_0_C1)=Proj0{interf_dir_parallel};  %first term of refinement formula in Lemma 2
   C(ind_start_int_dofs_patch_ref:ind_end_int_dofs_patch_ref,ndof_interior_C1+1:ndof_interior_C1+ndof_0_C1)=...
       Aux(setdiff(1:ndof_Bsp_ref,union(ind0_ref,ind1_ref)),ndof_interior_C1+1:ndof_interior_C1+ndof_0_C1);   %second term of refinement formula in Lemma 2
    
   %Taking care of the parts corresponding to Lemma 3 (Mario's notes on refinement)
   
   C(ndof_interior_C1_ref+ndof_0_C1_ref+1:end,ndof_interior_C1+ndof_0_C1+1:end)=Proj1{interf_dir_parallel};  %first term of refinement formula in Lemma 3
   C(ind_start_int_dofs_patch_ref:ind_end_int_dofs_patch_ref,ndof_interior_C1+ndof_0_C1+1:end)=...
       Aux(setdiff(1:ndof_Bsp_ref,union(ind0_ref,ind1_ref)),ndof_interior_C1+ndof_0_C1+1:end);   %second term of refinement formula in Lemma 3
    
  end
  
%What's exactly the meaning of the 3rd input? Specifying only some
%functions to be expressed in the new basis? Yes, but we'll take care of it later

% We have to change all this part following Mario's notes, and using Proj, Proj0 and Proj1
elseif (nargin == 3)
  nrows = hspace.space_of_level(lev).ndof; ncols = hspace.space_of_level(lev-1).ndof;
  size_alloc = numel (ind_coarse) * prod (hspace.space_of_level(lev).sp_patch{1}.degree + 1);
  rows = zeros (size_alloc, 1); cols = rows; vals = rows;
  ncounter = 0;
  for ipatch = 1:npatch
    nc_in_patch = ncounter+1;
    [ind_ptc, ~, local_indices] = intersect (ind_coarse, hspace.space_of_level(lev-1).gnum{ipatch});
    sub_coarse = cell (ndim, 1);
    [sub_coarse{:}] = ind2sub ([hspace.space_of_level(lev-1).sp_patch{ipatch}.ndof_dir, 1], local_indices);

    for ii = 1:numel(ind_ptc)
      Proj = hspace.Proj{lev-1, ipatch};
      Caux = 1;
      for idim = 1:ndim
        Caux = kron (Proj{idim}(:,sub_coarse{idim}(ii)), Caux);
      end
      
      [ir, ic, iv] = find (Caux);
%       rows(ncounter+(1:numel(ir))) = hspace.space_of_level(lev).gnum{ipatch}(ir);
%       cols(ncounter+(1:numel(ir))) = ic*ind_ptc(ii);
      rows(ncounter+(1:numel(ir))) = ir;
      cols(ncounter+(1:numel(ir))) = local_indices(ii);
      vals(ncounter+(1:numel(ir))) = iv;
      ncounter = ncounter + numel (ir);
    end

    rows(nc_in_patch:ncounter) = hspace.space_of_level(lev).gnum{ipatch}(rows(nc_in_patch:ncounter));
    cols(nc_in_patch:ncounter) = hspace.space_of_level(lev-1).gnum{ipatch}(cols(nc_in_patch:ncounter));
  end

  rows = rows(1:ncounter); cols = cols(1:ncounter); vals = vals(1:ncounter);
%   C = sparse (rows, cols, vals, hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof);
  if (~isempty (rows))
    C = accumarray ([rows,cols], vals.', [], @min, 0, true);
  else
    C = sparse (nrows, ncols);
  end
  if (size (C,1) < nrows || size(C,2) < ncols)
    C(nrows, ncols) = 0;
  end
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  C(indices,:) = 0;
end

end
