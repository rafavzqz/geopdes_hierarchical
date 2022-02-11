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
% Copyright (C) 2019--2022 Cesare Bracco, Rafael Vazquez
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

function C = matrix_basis_change__ (hspace, lev, ind_coarse, ind_fine)

% TO DO: ALLOW TO USE rows, cols, vals

Proj = hspace.Proj(lev-1,:);
Proj0 = hspace.Proj0(lev-1,:);
Proj1 = hspace.Proj1(lev-1,:);

if (nargin == 4)
  C = subdivision_matrix_two_levels_C1__ (hspace.space_of_level(lev-1), hspace.space_of_level(lev), Proj, Proj0, Proj1, ind_coarse, ind_fine);
elseif (nargin == 3)
  C = subdivision_matrix_two_levels_C1__ (hspace.space_of_level(lev-1), hspace.space_of_level(lev), Proj, Proj0, Proj1, ind_coarse);  
else
  C = subdivision_matrix_two_levels_C1__ (hspace.space_of_level(lev-1), hspace.space_of_level(lev), Proj, Proj0, Proj1);
end

if (hspace.truncated)
  indices = union (hspace.active{lev}, hspace.deactivated{lev});
  if (nargin == 4)
    [~,~,indices] = intersect (indices,ind_fine);
  end
  C(indices,:) = 0;
end

end

function C = subdivision_matrix_two_levels_C1__ (sp_coarse, sp_fine, Proj, Proj0, Proj1, varargin)

if (isa (sp_coarse.sp_patch{1}, 'sp_scalar'))
  is_scalar = true;
  ndim = size (Proj{1}, 2);
  all_dir = 1:ndim; % possible directions
% elseif (isa (sp_coarse.sp_patch{1}, 'sp_vector'))
%   is_scalar = false;
%   ndim = size (Proj{1}, 2);
else
  error ('Unknown space type')
end

% Number of patches, interfaces and vertices
% npatch = sp_coarse.npatch;
% nvert = numel (sp_coarse.vertices);

if (nargin <= 6)
  C = sparse (sp_fine.ndof, sp_coarse.ndof); %zero sparse matrix 
elseif (nargin == 7)
  C = sparse (numel(varargin{2}), numel(varargin{1}));
end

% Computation of the matrix for patch interior basis functions
  C = subdivision_interior (sp_coarse, sp_fine, Proj, varargin{:});

% Computation of the matrix for edge basis functions
  C = C + subdivision_edges (sp_coarse, sp_fine, Proj, Proj0, Proj1, varargin{:});
  
% Computation of the matrix for vertex basis functions
  C = C + subdivision_vertices (sp_coarse, sp_fine, Proj, Proj0, Proj1, varargin{:});

end

function C = subdivision_interior (sp_coarse, sp_fine, Proj, ind_coarse, ind_fine)

  npatch = sp_coarse.npatch;

  if (nargin == 5)
    C = sparse (numel(ind_fine), numel(ind_coarse));
    for iptc = 1:npatch
      spc_patch = sp_coarse.sp_patch{iptc};
      spf_patch = sp_fine.sp_patch{iptc};
    
% Compute the indices of interior B-splines
      interior_inds_coarse = sp_coarse.dofs_on_patch{iptc};
      interior_inds_fine = sp_fine.dofs_on_patch{iptc};

      [~,local_indices,ind_c] = intersect (interior_inds_coarse, ind_coarse);
      ind_coarse_on_patch = sp_coarse.interior_dofs_per_patch{iptc}(local_indices);
      [~,local_indices,ind_f] = intersect (interior_inds_fine, ind_fine);
      ind_fine_on_patch = sp_fine.interior_dofs_per_patch{iptc}(local_indices);
      
      Cpatch = subdivision_matrix_two_levels__ (spc_patch, spf_patch, Proj{iptc}, ind_coarse_on_patch, ind_fine_on_patch); 
      C(ind_f, ind_c) = Cpatch;
    end
  elseif (nargin == 4 || nargin == 3)
    C = sparse (sp_fine.ndof, sp_coarse.ndof);
    for iptc = 1:npatch
      spc_patch = sp_coarse.sp_patch{iptc};
      spf_patch = sp_fine.sp_patch{iptc};
      
% Compute the indices of interior B-splines
      ndof_dir_spn = spc_patch.ndof_dir;
      ndof_Bsp = prod(ndof_dir_spn);
      ndof_dir_spn_ref = spf_patch.ndof_dir; 
      ndof_Bsp_ref = prod(ndof_dir_spn_ref);
      boundary = struct(spc_patch.boundary);
      ind = union ([boundary.dofs], [boundary.adjacent_dofs]);
      boundary = struct(spf_patch.boundary);
      ind_ref = union ([boundary.dofs], [boundary.adjacent_dofs]);
      Bsp_indices_coarse = setdiff (1:ndof_Bsp, ind);
      Bsp_indices_fine = setdiff (1:ndof_Bsp_ref, ind_ref);
      interior_inds_coarse = sp_coarse.dofs_on_patch{iptc};
      interior_inds_fine = sp_fine.dofs_on_patch{iptc};

      if (nargin == 4)
        [~,local_indices,~] = intersect (interior_inds_coarse, ind_coarse);
        ind_coarse_on_patch = sp_coarse.interior_dofs_per_patch{iptc}(local_indices);
        Cpatch = subdivision_matrix_two_levels__ (spc_patch, spf_patch, Proj{iptc}, ind_coarse_on_patch);
      elseif (nargin == 3)
        Cpatch = subdivision_matrix_two_levels__ (spc_patch, spf_patch, Proj{iptc});
      end
      C(interior_inds_fine,interior_inds_coarse) = Cpatch(Bsp_indices_fine,Bsp_indices_coarse);
    end
  end
  
end

function C = subdivision_edges (sp_coarse, sp_fine, Proj, Proj0, Proj1, ind_coarse, ind_fine)

  if (nargin < 7)
    C = sparse (sp_fine.ndof, sp_coarse.ndof);
  elseif (nargin == 7)
    C = sparse (numel(ind_fine), numel(ind_coarse));
  end

  ndim = 2; 
  all_dir = 1:ndim;
  
  interfaces = sp_coarse.interfaces;
  nint = numel (interfaces);
  
  shift_inds_ref = cumsum ([0 cellfun(@numel, sp_fine.interior_dofs_per_patch)]);
  shift_inds_e = cumsum([0 cellfun(@numel, sp_coarse.dofs_on_edge)]);
  shift_inds_e_ref = cumsum([0 cellfun(@numel, sp_fine.dofs_on_edge)]);
  ndof_interior_C1 = sp_coarse.ndof_interior;
  ndof_interior_C1_ref = sp_fine.ndof_interior;
  for ii = 1:nint
    if (nargin == 6 || nargin == 7)
      all_inds_on_edge = sp_coarse.dofs_on_edge{ii};
      if (isempty (intersect (all_inds_on_edge, ind_coarse)))
        continue
      end
    end
    clear side_on_int 
    clear patch_on_int
    if (~isempty(interfaces(ii).patch1))
      side_on_int(1) = interfaces(ii).side1;
      patch_on_int(1) = interfaces(ii).patch1;
    end
    if (~isempty(interfaces(ii).patch2))
      side_on_int(2) = interfaces(ii).side2;
      patch_on_int(2) = interfaces(ii).patch2;
    end
    
    for ipatch = 1:length(patch_on_int) %local patch index
      side = side_on_int(ipatch);
      interf_dir_orthogonal = ceil (side/2);  %direction orthogonal to the interface
      interf_dir_parallel = setdiff (all_dir, interf_dir_orthogonal);  %direction(s) parallel to the interface
      interf_dir = interf_dir_parallel;
      degree = sp_coarse.sp_patch{patch_on_int(ipatch)}.degree(interf_dir);

      % Number of univariate B-splines on the interface: sp0 (p,r+1) and sp1 (p-1,r)
      ndof_0_C1 = length (sp_coarse.knots0_patches{patch_on_int(ipatch)}{interf_dir}) - degree - 1;
      ndof_0_C1_ref = length (sp_fine.knots0_patches{patch_on_int(ipatch)}{interf_dir}) - degree - 1;
      ndof_1_C1 = length (sp_coarse.knots1_patches{patch_on_int(ipatch)}{interf_dir}) - degree ;
      ndof_1_C1_ref = length (sp_fine.knots1_patches{patch_on_int(ipatch)}{interf_dir}) - degree;
        
      % Indices corresponding to the standard B-splines spanning edge functions
      boundary = struct(sp_fine.sp_patch{patch_on_int(ipatch)}.boundary);
      ind_ref = union ([boundary.dofs], [boundary.adjacent_dofs]);
        
      %define the 1D B-spline space parallel to the interface: spn (p,r)
      ndof_dir_spn_ref = sp_fine.sp_patch{patch_on_int(ipatch)}.ndof_dir;
      ndof_Bsp_ref = prod(ndof_dir_spn_ref);
        
      Bsp_indices_fine = setdiff (1:ndof_Bsp_ref, ind_ref); %Indices of interior functions in the finer level
      interior_inds_fine = shift_inds_ref(patch_on_int(ipatch))+1:shift_inds_ref(patch_on_int(ipatch)+1);

      % Auxiliary matrices (standard refinement matrix)
      Lambda = 1;
      Proj_patch = Proj{patch_on_int(ipatch)};
      for idim = 1:ndim
        Lambda = kron (Proj_patch{idim}, Lambda);
      end
      Proj0_patch = Proj0{patch_on_int(ipatch)};
      Proj1_patch = Proj1{patch_on_int(ipatch)};
        
      % Indices0 and indices1 respectively correspond to trace and derivative edge functions
      indices0_coarse = ndof_interior_C1 + shift_inds_e(ii) + (1:ndof_0_C1-6);
      indices0_fine = ndof_interior_C1_ref + shift_inds_e_ref(ii) + (1:ndof_0_C1_ref-6);
      indices1_coarse = ndof_interior_C1 + shift_inds_e(ii) + ndof_0_C1-6 + (1:ndof_1_C1-4);
      indices1_fine = ndof_interior_C1_ref + shift_inds_e_ref(ii) + ndof_0_C1_ref-6 + (1:ndof_1_C1_ref-4);
      
      if (nargin == 5)
        if (ipatch == 1) % Doing it this way, I don't need to care about the orientation
          C(indices0_fine,indices0_coarse) = Proj0_patch{interf_dir_parallel}(4:end-3,4:end-3);  %first term of refinement formula in Lemma 2
          C(indices1_fine,indices1_coarse) = (1/2)*Proj1_patch{interf_dir_parallel}(3:end-2,3:end-2);  %first term of refinement formula in Lemma 3
        end
        Aux = Lambda * sp_coarse.Cpatch{patch_on_int(ipatch)};
        C(interior_inds_fine,indices0_coarse) = Aux(Bsp_indices_fine,indices0_coarse);   %second term of refinement formula in Lemma 2
        C(interior_inds_fine,indices1_coarse) = Aux(Bsp_indices_fine,indices1_coarse);   %second term of refinement formula in Lemma 3

      elseif (nargin == 6)
        [ind_coarse0_on_edge,local_indices0,~] = intersect (indices0_coarse, ind_coarse);
        [ind_coarse1_on_edge,local_indices1,~] = intersect (indices1_coarse, ind_coarse);
        if (ipatch == 1) % Doing it this way, I don't need to care about the orientation
          C(indices0_fine,ind_coarse0_on_edge) = Proj0_patch{interf_dir_parallel}(4:end-3,local_indices0+3);  %first term of refinement formula in Lemma 2
          C(indices1_fine,ind_coarse1_on_edge) = (1/2)*Proj1_patch{interf_dir_parallel}(3:end-2,local_indices1+2);  %first term of refinement formula in Lemma 3 %WARNING: MULTIPLIED BY 1/2?
        end
        Aux = Lambda * sp_coarse.Cpatch{patch_on_int(ipatch)};
        C(interior_inds_fine,ind_coarse0_on_edge) = Aux(Bsp_indices_fine,ind_coarse0_on_edge);   %second term of refinement formula in Lemma 2
        C(interior_inds_fine,ind_coarse1_on_edge) = Aux(Bsp_indices_fine,ind_coarse1_on_edge);   %second term of refinement formula in Lemma 3

      elseif (nargin == 7)
        [ind_coarse0_on_edge,local_indices_c0,ind_c0] = intersect (indices0_coarse, ind_coarse);
        [ind_coarse1_on_edge,local_indices_c1,ind_c1] = intersect (indices1_coarse, ind_coarse);
        [~,local_indices_f0,ind_f0] = intersect (indices0_fine, ind_fine);
        [~,local_indices_f1,ind_f1] = intersect (indices1_fine, ind_fine);
        
        [~,local_indices_int,ind_int_f] = intersect (interior_inds_fine, ind_fine);
        if (ipatch == 1) % Doing it this way, I don't need to care about the orientation
          C(ind_f0,ind_c0) = Proj0_patch{interf_dir_parallel}(local_indices_f0+3,local_indices_c0+3);  %first term of refinement formula in Lemma 2
          C(ind_f1,ind_c1) = (1/2)*Proj1_patch{interf_dir_parallel}(local_indices_f1+2,local_indices_c1+2);  %first term of refinement formula in Lemma 3 %WARNING: MULTIPLIED BY 1/2?
        end
        Aux = Lambda * sp_coarse.Cpatch{patch_on_int(ipatch)};
        C(ind_int_f,ind_c0) = Aux(Bsp_indices_fine(local_indices_int),ind_coarse0_on_edge);   %second term of refinement formula in Lemma 2
        C(ind_int_f,ind_c1) = Aux(Bsp_indices_fine(local_indices_int),ind_coarse1_on_edge);   %second term of refinement formula in Lemma 3
      end
    end
  end

end

function C = subdivision_vertices (sp_coarse, sp_fine, Proj, Proj0, Proj1, ind_coarse, ind_fine)

  if (nargin < 7)
    C = sparse (sp_fine.ndof, sp_coarse.ndof);
  elseif (nargin == 7)
    C = sparse (numel(ind_fine), numel(ind_coarse));
  end

  ndim = 2;
  nvert = numel (sp_coarse.vertices);
  
% Auxiliary arrays to easily compute the direction of the interface,
% according to patch_reorientation, without using three nested if-else.
  interf_dirs_patch1(:,:,1) = [1 2; 2 1];
  interf_dirs_patch1(:,:,2) = [2 1; 1 2];
  interf_dirs_patch2 = interf_dirs_patch1(:,:,[2 1]);

  for iv = 1:nvert
    %global indices of functions associated to this vertex
    indices_v_coarse = sp_coarse.dofs_on_vertex{iv};
    indices_v_fine = sp_fine.dofs_on_vertex{iv};
    
    if ((nargin == 6 || nargin == 7) && isempty (intersect(indices_v_coarse, ind_coarse)))
      continue
    end
    patches = sp_coarse.vertices(iv).patches;
    edges = sp_coarse.vertices(iv).edges;
    
    %a) part of the matrix describing the dependence on the same vertex functions on the finer level
    sigma_coarse = sp_coarse.vertex_function_matrices{1,iv};
    sigma_fine = sp_fine.vertex_function_matrices{1,iv};
    sigma_rat = sigma_coarse/sigma_fine;
    sigma_vec = [1 sigma_rat sigma_rat^2 sigma_rat sigma_rat^2 sigma_rat^2];
    Dsigma = diag(sigma_vec);
    if (nargin == 5)
      C(indices_v_fine,indices_v_coarse) = Dsigma;
    elseif (nargin == 6)
      [ind_coarse_on_vertex,local_indices,~] = intersect (indices_v_coarse, ind_coarse);
      C(indices_v_fine, ind_coarse_on_vertex) = Dsigma(:,local_indices);
    elseif (nargin == 7)
      [~,local_indices_c,ind_c] = intersect (indices_v_coarse, ind_coarse);
      [~,local_indices_f,ind_f] = intersect (indices_v_fine, ind_fine);        
      C(ind_f, ind_c) = Dsigma(local_indices_f,local_indices_c);
    end
    
    %b) part of the matrix describing the dependence on edge and interior functions on the finer level
    for ip = 1:numel(patches)
      ip_plus_1 = mod(ip, sp_coarse.vertices(iv).valence_e) + 1;
        
      %determining interf_dir (interface direction)
      aux_indices = num2cell (sp_coarse.vertices(iv).patch_reorientation(ip,:) + 1);
      interf_dir = interf_dirs_patch1(aux_indices{:});
        
      % Get the indices of the interior standard B-splines (finer level)
      ndof_dir_spn_ref = sp_fine.sp_patch{patches(ip)}.ndof_dir;
      [indx, indy] = ndgrid (3:ndof_dir_spn_ref(1)-2, 3:ndof_dir_spn_ref(2)-2);
      int_ref = sub2ind (ndof_dir_spn_ref, indx(:).', indy(:).');
        
      K_prev = sp_coarse.vertex_function_matrices{2,iv}{ip}.K_prev;
      E_prev = sp_coarse.vertex_function_matrices{2,iv}{ip}.E_prev;
      if (sp_coarse.vertices(iv).edge_orientation(ip) == -1) 
        E_prev(:,[4 5]) = -E_prev(:,[4 5]);
        E_prev = E_prev(:,[3 2 1 5 4]);
      end
      K_next = sp_coarse.vertex_function_matrices{2,iv}{ip}.K_next;
      E_next = sp_coarse.vertex_function_matrices{2,iv}{ip}.E_next;
      if (sp_coarse.vertices(iv).edge_orientation(ip_plus_1) == -1) 
        E_next(:,[4 5]) = -E_next(:,[4 5]);
        E_next = E_next(:,[3 2 1 5 4]);
        K_next([4 5],:) = -K_next([4 5],:);
        K_next = K_next([3 2 1 5 4],:);
      end
        
      %Auxiliary matrices (standard refinement matrix)
      Lambda = 1;
      Proj_patch = Proj{patches(ip)};
      for idim = 1:ndim
        Lambda = kron (Proj_patch{idim}, Lambda);
      end
      Proj0_patch = Proj0{patches(ip)};
      Proj1_patch = Proj1{patches(ip)};
      dim_sp0 = size (Proj0_patch{interf_dir},2);
      dim_sp1 = size (Proj1_patch{interf_dir},2);
      dim_sp0_ref = size (Proj0_patch{interf_dir},1);
      dim_sp1_ref = size (Proj1_patch{interf_dir},1);

      % Dependence on edge functions (using previous edge)
      inactive_edge = [1 2 3 dim_sp0+1 dim_sp0+2];
      active_edge_ref = [4:dim_sp0_ref-3 dim_sp0_ref+3:dim_sp0_ref+dim_sp1_ref-2];
      if (sp_coarse.vertices(iv).edge_orientation(ip) == 1)
        Aux_edge_disc = [Proj0_patch{interf_dir},          zeros(dim_sp0_ref,dim_sp1);...
                         zeros(dim_sp1_ref,dim_sp0), (1/2)*Proj1_patch{interf_dir}];
        ind_edge_ref = sp_fine.dofs_on_edge{edges(ip)};
      else
        Aux_edge_disc = [Proj0_patch{interf_dir},          zeros(dim_sp0_ref,dim_sp1);...
                       zeros(dim_sp1_ref,dim_sp0), (-1/2)*Proj1_patch{interf_dir}];
        nn0 = numel(4:dim_sp0_ref-3);
        nn1 = numel(3:dim_sp1_ref-2);
        ind_edge_ref = sp_fine.dofs_on_edge{edges(ip)}([nn0:-1:1, nn0+(nn1:-1:1)]);
      end

      if (nargin == 5)
        C(ind_edge_ref,indices_v_coarse) = Aux_edge_disc(active_edge_ref,inactive_edge)*K_prev;
      elseif (nargin == 6)
        [ind_coarse_on_vertex,local_indices,~] = intersect (indices_v_coarse, ind_coarse);
        C(ind_edge_ref,ind_coarse_on_vertex) = Aux_edge_disc(active_edge_ref,inactive_edge)*K_prev(:,local_indices);
      elseif (nargin == 7)
        [~,local_indices_c,ind_c] = intersect (indices_v_coarse, ind_coarse);
        [~,local_indices_f,ind_f] = intersect (ind_edge_ref, ind_fine);
        C(ind_f, ind_c) = Aux_edge_disc(active_edge_ref(local_indices_f),inactive_edge)*K_prev(:,local_indices_c);
      end
        
      % Dependence on standard B-splines, first through edge functions (K matrices) and then directly (V matrix)
      if (sp_coarse.vertices(iv).edge_orientation(ip) == -1)
        K_prev([4 5],:) = -K_prev([4 5],:);
        K_prev = K_prev([3 2 1 5 4],:);
      end
      Aux_prev = Lambda * E_prev;
      Aux_next = Lambda * E_next;
      V = sp_coarse.vertex_function_matrices{2,iv}{ip}.V;
      Aux = Lambda*V;
      ind_int_ref = sp_fine.dofs_on_patch{patches(ip)}; %shift_inds_ref(patches(ip))+1:shift_inds_ref(patches(ip)+1);
      
      if (nargin == 5)
        C(ind_int_ref,indices_v_coarse) = Aux_prev(int_ref,:)*K_prev + Aux_next(int_ref,:)*K_next;
        C(ind_int_ref,indices_v_coarse) = C(ind_int_ref,indices_v_coarse) - Aux(int_ref,:);
      elseif (nargin == 6)
        [ind_coarse_on_vertex,local_indices,~] = intersect (indices_v_coarse, ind_coarse);
        C(ind_int_ref,ind_coarse_on_vertex) = Aux_prev(int_ref,:)*K_prev(:,local_indices) + Aux_next(int_ref,:)*K_next(:,local_indices);
        C(ind_int_ref,indices_v_coarse) = C(ind_int_ref,indices_v_coarse) - Aux(int_ref,local_indices);
      elseif (nargin == 7)
        [~,local_indices_c,ind_c] = intersect (indices_v_coarse, ind_coarse);
        [~,local_indices_f,ind_f] = intersect (ind_int_ref, ind_fine);
        C(ind_f,ind_c) = Aux_prev(int_ref(local_indices_f),:)*K_prev(:,local_indices_c) + Aux_next(int_ref(local_indices_f),:)*K_next(:,local_indices_c);
        C(ind_f,ind_c) = C(ind_f,ind_c) - Aux(int_ref(local_indices_f),local_indices_c);
      end
    end
    
   %if there are MORE EDGES THAN PATCHES (boundary vertex), we must find
   %the dependence on the corresponding edge functions
    if (numel(edges) > numel(patches))
      ip = numel(patches); %we are on the last patch and we consider the (ip+1)-th edge
      %determining interf_dir
      aux_indices = num2cell (sp_coarse.vertices(iv).patch_reorientation(ip,:) + 1);
      interf_dir = interf_dirs_patch2(aux_indices{:});
        
      K_next=sp_coarse.vertex_function_matrices{2,iv}{ip}.K_next;
        
      %Auxiliary matrices (standard refinement matrix)
      Lambda = 1;
      Proj_patch = Proj{patches(ip)};
      for idim = 1:ndim
        Lambda = kron (Proj_patch{idim}, Lambda);
      end
      Proj0_patch = Proj0{patches(ip)};
      Proj1_patch = Proj1{patches(ip)};

      dim_sp0=size(Proj0_patch{interf_dir},2);
      dim_sp1=size(Proj1_patch{interf_dir},2);
      dim_sp0_ref=size(Proj0_patch{interf_dir},1);
      dim_sp1_ref=size(Proj1_patch{interf_dir},1);

      % Dependence on edge functions (using next edge)
      inactive_edge = [1 2 3 dim_sp0+1 dim_sp0+2];
      active_edge_ref = [4:dim_sp0_ref-3 dim_sp0_ref+3:dim_sp0_ref+dim_sp1_ref-2];
      if (sp_coarse.vertices(iv).edge_orientation(ip+1)==1)
        Aux_edge_disc = [Proj0_patch{interf_dir},          zeros(dim_sp0_ref,dim_sp1);...
                         zeros(dim_sp1_ref,dim_sp0), (1/2)*Proj1_patch{interf_dir}];
        ind_edge_ref = sp_fine.dofs_on_edge{edges(ip+1)};
      else
        Aux_edge_disc = [Proj0_patch{interf_dir},          zeros(dim_sp0_ref,dim_sp1);...
                         zeros(dim_sp1_ref,dim_sp0), (-1/2)*Proj1_patch{interf_dir}];
        nn0 = numel(4:dim_sp0_ref-3);
        nn1 = numel(3:dim_sp1_ref-2);
        ind_edge_ref = sp_fine.dofs_on_edge{edges(ip+1)}([nn0:-1:1, nn0+(nn1:-1:1)]);
      end
      if (nargin == 5)
        C(ind_edge_ref,indices_v_coarse) = Aux_edge_disc(active_edge_ref,inactive_edge)*K_next;
      elseif (nargin == 6)
        [ind_coarse_on_vertex,local_indices,~] = intersect (indices_v_coarse, ind_coarse);
        C(ind_edge_ref,ind_coarse_on_vertex) = Aux_edge_disc(active_edge_ref,inactive_edge)*K_next(:,local_indices);
      elseif (nargin == 7)
        [~,local_indices_c,ind_c] = intersect (indices_v_coarse, ind_coarse);
        [~,local_indices_f,ind_f] = intersect (ind_edge_ref, ind_fine);
        C(ind_f, ind_c) = Aux_edge_disc(active_edge_ref(local_indices_f),inactive_edge)*K_next(:,local_indices_c);
      end
    end 
  end

end