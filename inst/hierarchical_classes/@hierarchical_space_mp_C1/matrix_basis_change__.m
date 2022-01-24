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

function C = matrix_basis_change__ (hspace, lev)




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
  all_dir = 1:ndim; % possible directions
% elseif (isa (hspace.space_of_level(1).sp_patch{1}, 'sp_vector'))
%   is_scalar = false;
%   ndim = size (hspace.Proj{1}, 2);
else
  error ('Unknown space type')
end

%MP update:
npatch = hspace.space_of_level(1).npatch;  %number of patches
nint = numel(hspace.space_of_level(1).interfaces);  %number of interfaces
nvert = numel(hspace.space_of_level(1).vertices);  %number of vertices

% For now we only use B-splines, and scalar spaces

% We have to change all this part following Mario's notes, and using Proj, Proj0 and Proj1
% How are the three types of basis functions organized (indexed) in the hspace object?
if (nargin < 3)
    
  C = sparse (hspace.space_of_level(lev).ndof, hspace.space_of_level(lev-1).ndof); %zero sparse matrix 

  %Cesare's note: this defines the dimension of the matrix C (considers full multipatch spaces on consecutive levels)
  % in the C^0 case, for each patch you compute the classical matrix going
  % from level l-1 to level l (tensor product basis), and then put it in
  % the right place in the big matrix, according to the global numbering of
  % the dofs of the patch.
  % What does Proj (and Proj0, Proj1) contain in our case?

% Shifting indices for interior degrees of freedom
  shift_inds = cumsum ([0 cellfun(@numel, hspace.space_of_level(lev-1).interior_dofs_per_patch)]);
  shift_inds_ref = cumsum ([0 cellfun(@numel, hspace.space_of_level(lev).interior_dofs_per_patch)]);
  shift_inds_e = cumsum([0 cellfun(@numel, hspace.space_of_level(lev-1).dofs_on_edge)]);
  shift_inds_e_ref = cumsum([0 cellfun(@numel, hspace.space_of_level(lev).dofs_on_edge)]);
  shift_inds_v = cumsum([0 cellfun(@numel, hspace.space_of_level(lev-1).dofs_on_vertex)]);
  shift_inds_v_ref = cumsum([0 cellfun(@numel, hspace.space_of_level(lev).dofs_on_vertex)]);
  
%in Cpatch{iptc}, columns 1:sp.ndof_interior represent the regular (interior) B-splines (all patches),  
  ndof_interior_C1 = hspace.space_of_level(lev-1).ndof_interior;
  ndof_interior_C1_ref = hspace.space_of_level(lev).ndof_interior; %same for the finer level
  ndof_edge_C1 = hspace.space_of_level(lev-1).ndof_edges;
  ndof_edge_C1_ref = hspace.space_of_level(lev).ndof_edges;
%   ndof_vertices_C1 = hspace.space_of_level(lev-1).ndof_vertices;
%   ndof_vertices_C1_ref = hspace.space_of_level(lev).ndof_vertices;
  Cpatch = hspace.space_of_level(lev-1).Cpatch;  
  
  for ip=1:npatch
    Lambda = 1;
    Proj = hspace.Proj{lev-1, ip};
    for idim = 1:ndim
      Lambda = kron (Proj{idim}, Lambda);
    end

    %Dimension of standard bivariate spline spaces on the patch
    ndof_dir_spn = hspace.space_of_level(lev-1).sp_patch{ip}.ndof_dir;
    ndof_Bsp = prod(ndof_dir_spn);
    ndof_dir_spn_ref = hspace.space_of_level(lev).sp_patch{ip}.ndof_dir; 
    ndof_Bsp_ref = prod(ndof_dir_spn_ref);

    %Indices of the B-splines which are not active (not interior)
    boundary = struct(hspace.space_of_level(lev-1).sp_patch{ip}.boundary);
    ind = union ([boundary.dofs], [boundary.adjacent_dofs]);
    boundary = struct(hspace.space_of_level(lev).sp_patch{ip}.boundary);
    ind_ref = union ([boundary.dofs], [boundary.adjacent_dofs]);
%     ind=sub2ind(ndof_dir_spn,[1*ones(1,ndof_dir_spn(2)) 2*ones(1,ndof_dir_spn(2))...
%         (ndof_dir_spn(1)-1)*ones(1,ndof_dir_spn(2)) ndof_dir_spn(1)*ones(1,ndof_dir_spn(2))],repmat(1:ndof_dir_spn(2),1,4));
%     ind=[ind sub2ind(ndof_dir_spn,repmat(2:(ndof_dir_spn(1)-1),1,4),[1*ones(1,ndof_dir_spn(1)-2) 2*ones(1,ndof_dir_spn(1)-2)...
%         (ndof_dir_spn(2)-1)*ones(1,ndof_dir_spn(1)-2) ndof_dir_spn(2)*ones(1,ndof_dir_spn(1)-2)])];
%     ind_ref=sub2ind(ndof_dir_spn_ref,[1*ones(1,ndof_dir_spn_ref(2)) 2*ones(1,ndof_dir_spn_ref(2))...
%         (ndof_dir_spn_ref(1)-1)*ones(1,ndof_dir_spn_ref(2)) ndof_dir_spn_ref(1)*ones(1,ndof_dir_spn_ref(2))],repmat(1:ndof_dir_spn_ref(2),1,4));
%     ind_ref=[ind_ref sub2ind(ndof_dir_spn_ref,repmat(2:(ndof_dir_spn_ref(1)-1),1,4),[1*ones(1,ndof_dir_spn_ref(1)-2) 2*ones(1,ndof_dir_spn_ref(1)-2)...
%         (ndof_dir_spn_ref(2)-1)*ones(1,ndof_dir_spn_ref(1)-2) ndof_dir_spn_ref(2)*ones(1,ndof_dir_spn_ref(1)-2)])];
    %Refinement of interior functions
    Bsp_indices_coarse = setdiff (1:ndof_Bsp, ind);
    Bsp_indices_fine = setdiff (1:ndof_Bsp_ref, ind_ref);
    interior_inds_coarse = shift_inds(ip)+1:shift_inds(ip+1);
    interior_inds_fine = shift_inds_ref(ip)+1:shift_inds_ref(ip+1);
    
    C(interior_inds_fine,interior_inds_coarse) = Lambda(Bsp_indices_fine, Bsp_indices_coarse);  %not all Lambda
  end
  
  interfaces = hspace.space_of_level(lev-1).interfaces;
  for ii = 1:nint
    %get interf_dir=1,2 according to the interface being horizontal or vertical
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
        degree = hspace.space_of_level(lev-1).sp_patch{patch_on_int(ipatch)}.degree(interf_dir);
        
        %define the 1D B-spline space parallel to the interface: spn (p,r)
        ndof_dir_spn = hspace.space_of_level(lev-1).sp_patch{patch_on_int(ipatch)}.ndof_dir;
        ndof_spn = ndof_dir_spn(interf_dir);
        ndof_Bsp = prod(ndof_dir_spn);  %dimension of bivariate space
        ndof_dir_spn_ref = hspace.space_of_level(lev).sp_patch{patch_on_int(ipatch)}.ndof_dir; %same for the finer level
        ndof_spn_ref = ndof_dir_spn_ref(interf_dir);     %same for the finer level
        ndof_Bsp_ref = prod(ndof_dir_spn_ref);   %same for the finer level
        
        % Get the indices ind0, ind1 corresponding to the standard B-splines spanning edge functions
        boundary = struct(hspace.space_of_level(lev-1).sp_patch{patch_on_int(ipatch)}.boundary);
        ind = union ([boundary.dofs], [boundary.adjacent_dofs]);
        boundary = struct(hspace.space_of_level(lev).sp_patch{patch_on_int(ipatch)}.boundary);
        ind_ref = union ([boundary.dofs], [boundary.adjacent_dofs]);
%         ind=sub2ind(ndof_dir_spn,[1*ones(1,ndof_dir_spn(2)) 2*ones(1,ndof_dir_spn(2))...
%         (ndof_dir_spn(1)-1)*ones(1,ndof_dir_spn(2)) ndof_dir_spn(1)*ones(1,ndof_dir_spn(2))],repmat(1:ndof_dir_spn(2),1,4));
%         ind=[ind sub2ind(ndof_dir_spn,repmat(2:(ndof_dir_spn(1)-1),1,4),[1*ones(1,ndof_dir_spn(1)-2) 2*ones(1,ndof_dir_spn(1)-2)...
%         (ndof_dir_spn(2)-1)*ones(1,ndof_dir_spn(1)-2) ndof_dir_spn(2)*ones(1,ndof_dir_spn(1)-2)])];  
%         ind_ref=sub2ind(ndof_dir_spn_ref,[1*ones(1,ndof_dir_spn_ref(2)) 2*ones(1,ndof_dir_spn_ref(2))...
%         (ndof_dir_spn_ref(1)-1)*ones(1,ndof_dir_spn_ref(2)) ndof_dir_spn_ref(1)*ones(1,ndof_dir_spn_ref(2))],repmat(1:ndof_dir_spn_ref(2),1,4));
%         ind_ref=[ind_ref sub2ind(ndof_dir_spn_ref,repmat(2:(ndof_dir_spn_ref(1)-1),1,4),[1*ones(1,ndof_dir_spn_ref(1)-2) 2*ones(1,ndof_dir_spn_ref(1)-2)...
%         (ndof_dir_spn_ref(2)-1)*ones(1,ndof_dir_spn_ref(1)-2) ndof_dir_spn_ref(2)*ones(1,ndof_dir_spn_ref(1)-2)])];

        %define the 1D B-spline spaces parallel to the interface: sp0 (p,r+1) and sp1 (p-1,r)
        %(actually it is enough to know the number of degrees of freedom of these two spaces)
        ndof_0_C1 = length (hspace.space_of_level(lev-1).knots0_patches{patch_on_int(ipatch)}{interf_dir}) - degree - 1;
        ndof_0_C1_ref = length (hspace.space_of_level(lev).knots0_patches{patch_on_int(ipatch)}{interf_dir}) - degree - 1; %same for the finer level
        ndof_1_C1 = length (hspace.space_of_level(lev-1).knots1_patches{patch_on_int(ipatch)}{interf_dir}) - degree ;
        ndof_1_C1_ref = length (hspace.space_of_level(lev).knots1_patches{patch_on_int(ipatch)}{interf_dir}) - degree; %same for the finer level
        
        %Auxiliary matrices (standard refinement matrix)
        Lambda = 1;
        Proj = hspace.Proj{lev-1, patch_on_int(ipatch)};
        for idim = 1:ndim
          Lambda = kron (Proj{idim}, Lambda);
        end
        Proj0 = hspace.Proj0{lev-1, patch_on_int(ipatch)};
        Proj1 = hspace.Proj1{lev-1, patch_on_int(ipatch)};
        Aux = Lambda * Cpatch{patch_on_int(ipatch)};
        
        %Bsp_indices_coarse = setdiff (1:ndof_Bsp, ind);
        Bsp_indices_fine = setdiff (1:ndof_Bsp_ref, ind_ref); %Indices of interior functions in the finer levels
        %interior_inds_coarse = shift_inds(ipatch)+1:shift_inds(ipatch+1);
        interior_inds_fine = shift_inds_ref(patch_on_int(ipatch))+1:shift_inds_ref(patch_on_int(ipatch)+1);
        
        %Taking care of the parts corresponding to trace edge functions
        indices0_coarse = ndof_interior_C1 + shift_inds_e(ii) + (1:ndof_0_C1-6); %also to be modified (taking into account trace functions of previous patches)
        indices0_fine = ndof_interior_C1_ref + shift_inds_e_ref(ii) + (1:ndof_0_C1_ref-6);
        if (ipatch == 1) % Doing it this way, I don't need to care about the orientation
          C(indices0_fine,indices0_coarse) = Proj0{interf_dir_parallel}(4:end-3,4:end-3);  %first term of refinement formula in Lemma 2
        end
        C(interior_inds_fine,indices0_coarse) = Aux(Bsp_indices_fine,indices0_coarse);   %second term of refinement formula in Lemma 2

       %Taking care of the parts corresponding to derivative edge functions
        indices1_coarse = ndof_interior_C1 + shift_inds_e(ii) + ndof_0_C1-6 + (1:ndof_1_C1-4); %also to be modified (taking into account trace functions of previous patches)
        indices1_fine = ndof_interior_C1_ref + shift_inds_e_ref(ii) + ndof_0_C1_ref-6 + (1:ndof_1_C1_ref-4);
        if (ipatch == 1) % Doing it this way, I don't need to care about the orientation
          C(indices1_fine,indices1_coarse) = (1/2)*Proj1{interf_dir_parallel}(3:end-2,3:end-2);  %first term of refinement formula in Lemma 3 %WARNING: MULTIPLIED BY 1/2?
        end
        C(interior_inds_fine,indices1_coarse) = Aux(Bsp_indices_fine,indices1_coarse);   %second term of refinement formula in Lemma 3
    end
  end
  
  %For this part, do we need to save the sigmas, Ks and Vs in the
  %hierarchical space structure? YES
  for iv = 1:nvert
    
    %refinement of vertex functions
    patches = hspace.space_of_level(lev-1).vertices(iv).patches;
    edges = hspace.space_of_level(lev-1).vertices(iv).edges;
    
    %global indices of this vertex functions (coarse level)
%     indices_v_coarse = hspace.space_of_level(lev-1).dofs_on_vertex{iv}
    indices_v_coarse = ndof_interior_C1 + ndof_edge_C1 + shift_inds_v(iv) + (1:6);
    indices_v_ref = ndof_interior_C1_ref + ndof_edge_C1_ref + shift_inds_v_ref(iv) + (1:6);
    
    %a) part of the matrix describing the dependence on the same vertex functions on the finer level
    sigma_coarse = hspace.space_of_level(lev-1).vertex_function_matrices{1,iv};
    sigma_fine = hspace.space_of_level(lev).vertex_function_matrices{1,iv};
    sigma_rat = sigma_coarse/sigma_fine;
    sigma_vec = [1 sigma_rat sigma_rat^2 sigma_rat sigma_rat^2 sigma_rat^2];
    C(indices_v_ref,indices_v_coarse) = diag(sigma_vec);
    
    %b) part of the matrix describing the dependence on edge and interior functions on the finer level
    for ip = 1:numel(patches)
        ip_plus_1=mod(ip, hspace.space_of_level(lev-1).vertices(iv).valence_e) + 1;
        prev_edge = edges(ip); %global index of the previous edge
        next_edge = edges(mod(ip, hspace.space_of_level(lev-1).vertices(iv).valence_e) + 1); %global index of the next edge (wrong, but not used)
        
        %determining interf_dir (interface direction)
        if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,3)==0 %x and y not inverted between them
            if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,1)==0   
                if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,2)==0 
                    interf_dir=1;
                else
                    interf_dir=2;
                end
            else
                if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,2)==0 
                    interf_dir=2;
                else
                    interf_dir=1;
                end
            end
        else
            if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,1)==0   
                if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,2)==0 
                    interf_dir=2;
                else
                    interf_dir=1;
                end
            else
                if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,2)==0 
                    interf_dir=1;
                else
                    interf_dir=2;
                end
            end            
        end
        
        % Get the indices of the interior standard B-splines (finer level)
        ndof_dir_spn_ref = hspace.space_of_level(lev).sp_patch{patches(ip)}.ndof_dir;
        [indx, indy] = ndgrid (3:ndof_dir_spn_ref(1)-2, 3:ndof_dir_spn_ref(2)-2);
        int_ref = sub2ind (ndof_dir_spn_ref, indx(:).', indy(:).');
        
        K_prev = hspace.space_of_level(lev-1).vertex_function_matrices{2,iv}{ip}.K_prev;
        E_prev = hspace.space_of_level(lev-1).vertex_function_matrices{2,iv}{ip}.E_prev;
        if hspace.space_of_level(lev-1).vertices(iv).edge_orientation(ip)==-1 
            E_prev(:,[4 5]) = -E_prev(:,[4 5]);
            E_prev = E_prev(:,[3 2 1 5 4]);
        end
        K_next = hspace.space_of_level(lev-1).vertex_function_matrices{2,iv}{ip}.K_next;
        E_next = hspace.space_of_level(lev-1).vertex_function_matrices{2,iv}{ip}.E_next;
        if hspace.space_of_level(lev-1).vertices(iv).edge_orientation(ip_plus_1)==-1 
            E_next(:,[4 5]) = -E_next(:,[4 5]);
            E_next = E_next(:,[3 2 1 5 4]);
            K_next([4 5],:) = -K_next([4 5],:);
            K_next = K_next([3 2 1 5 4],:);
        end
        
        %Auxiliary matrices (standard refinement matrix)
        Lambda = 1;
        Proj = hspace.Proj{lev-1, patches(ip)};
        for idim = 1:ndim
            Lambda = kron (Proj{idim}, Lambda);
        end
        Proj0 = hspace.Proj0{lev-1, patches(ip)};
        Proj1 = hspace.Proj1{lev-1, patches(ip)};
        dim_sp0 = size (Proj0{interf_dir},2);
        dim_sp1 = size (Proj1{interf_dir},2);
        dim_sp0_ref = size (Proj0{interf_dir},1);
        dim_sp1_ref = size (Proj1{interf_dir},1);
        
        %dependence on edge functions (finer level) of prev edge
        if (hspace.space_of_level(lev-1).vertices(iv).edge_orientation(ip)==1)
            inactive_edge = [1 2 3 dim_sp0+1 dim_sp0+2];
            active_edge_ref = [4:dim_sp0_ref-3 dim_sp0_ref+3:dim_sp0_ref+dim_sp1_ref-2];
            Aux_edge_disc = [Proj0{interf_dir},          zeros(dim_sp0_ref,dim_sp1);...
                             zeros(dim_sp1_ref,dim_sp0), (1/2)*Proj1{interf_dir}];
            ind_edge_ref = ndof_interior_C1_ref + ...
                           (shift_inds_e_ref(edges(ip))+1:shift_inds_e_ref(edges(ip)+1));
        else
            inactive_edge = [1 2 3 dim_sp0+1 dim_sp0+2];
            active_edge_ref = [4:dim_sp0_ref-3 dim_sp0_ref+3:dim_sp0_ref+dim_sp1_ref-2];
            Aux_edge_disc = [Proj0{interf_dir},          zeros(dim_sp0_ref,dim_sp1);...
                             zeros(dim_sp1_ref,dim_sp0), (-1/2)*Proj1{interf_dir}];
            nn0 = numel(4:dim_sp0_ref-3);
            nn1 = numel(3:dim_sp1_ref-2);
            ind_edge_ref = ndof_interior_C1_ref + shift_inds_e_ref(edges(ip)) + ...
                            [nn0:-1:1, nn0+(nn1:-1:1)];
        end
        C(ind_edge_ref,indices_v_coarse) = Aux_edge_disc(active_edge_ref,inactive_edge)*K_prev;
        
        %dependence on standard B-splines (finer level) through coarse edge functions of prev and next edge
        if hspace.space_of_level(lev-1).vertices(iv).edge_orientation(ip)==-1 
            K_prev([4 5],:) = -K_prev([4 5],:);
            K_prev = K_prev([3 2 1 5 4],:);
        end
        Aux_prev = Lambda * E_prev;
        Aux_next = Lambda * E_next;
        ind_int_ref = shift_inds_ref(patches(ip))+1:shift_inds_ref(patches(ip)+1);
        C(ind_int_ref,indices_v_coarse) = Aux_prev(int_ref,:)*K_prev+Aux_next(int_ref,:)*K_next;
        
        %dependence on standard B-splines (finer level)
        V = hspace.space_of_level(lev-1).vertex_function_matrices{2,iv}{ip}.V;
        Aux = Lambda*V;
        C(ind_int_ref,indices_v_coarse) = C(ind_int_ref,indices_v_coarse)-Aux(int_ref,:);
    end
    
   %if there are MORE EDGES THAN PATCHES (boundary vertex), we must find
   %the dependence on the corresponding edge functions
   if numel(edges)>numel(patches)
       ip=numel(patches); %we are on the last patch and we consider the (ip+1)-th edge
        %determining interf_dir
        if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,3)==0 %x and y not inverted between them
            if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,1)==0   
                if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,2)==0 
                    interf_dir=2;
                else
                    interf_dir=1;
                end
            else
                if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,2)==0 
                    interf_dir=1;
                else
                    interf_dir=2;
                end
            end
        else
            if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,1)==0   
                if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,2)==0 
                    interf_dir=1;
                else
                    interf_dir=2;
                end
            else
                if hspace.space_of_level(lev-1).vertices(iv).patch_reorientation(ip,2)==0 
                    interf_dir=2;
                else
                    interf_dir=1;
                end
            end            
        end
        
        K_next=hspace.space_of_level(lev-1).vertex_function_matrices{2,iv}{ip}.K_next;
        
        %Auxiliary matrices (standard refinement matrix)
        Lambda = 1;
        Proj = hspace.Proj{lev-1, patches(ip)};
        for idim = 1:ndim
            Lambda = kron (Proj{idim}, Lambda);
        end
        Proj0 = hspace.Proj0{lev-1, patches(ip)};
        Proj1 = hspace.Proj1{lev-1, patches(ip)};
        
        dim_sp0=size(Proj0{interf_dir},2);
        dim_sp1=size(Proj1{interf_dir},2);
        dim_sp0_ref=size(Proj0{interf_dir},1);
        dim_sp1_ref=size(Proj1{interf_dir},1);
        
        %dependence on edge functions (finer level) of next edge
        if (hspace.space_of_level(lev-1).vertices(iv).edge_orientation(ip+1)==1)
            inactive_edge = [1 2 3 dim_sp0+1 dim_sp0+2];
            active_edge_ref = [4:dim_sp0_ref-3 dim_sp0_ref+3:dim_sp0_ref+dim_sp1_ref-2];
            Aux_edge_disc = [Proj0{interf_dir},          zeros(dim_sp0_ref,dim_sp1);...
                             zeros(dim_sp1_ref,dim_sp0), (1/2)*Proj1{interf_dir}];
            ind_edge_ref = ndof_interior_C1_ref + ...
                           (shift_inds_e_ref(edges(ip+1))+1:shift_inds_e_ref(edges(ip+1)+1));
        else
            inactive_edge = [1 2 3 dim_sp0+1 dim_sp0+2];
            active_edge_ref = [4:dim_sp0_ref-3 dim_sp0_ref+3:dim_sp0_ref+dim_sp1_ref-2];
            Aux_edge_disc = [Proj0{interf_dir},          zeros(dim_sp0_ref,dim_sp1);...
                             zeros(dim_sp1_ref,dim_sp0), (-1/2)*Proj1{interf_dir}];
            nn0 = numel(4:dim_sp0_ref-3);
            nn1 = numel(3:dim_sp1_ref-2);
            ind_edge_ref = ndof_interior_C1_ref + shift_inds_e_ref(edges(ip+1)) + ...
                            [nn0:-1:1, nn0+(nn1:-1:1)];
        end       
        C(ind_edge_ref,indices_v_coarse)=Aux_edge_disc(active_edge_ref,inactive_edge)*K_next;
   end 
   
  end
  


  
% We have to change all this part following Mario's notes, and using Proj, Proj0 and Proj1
elseif (nargin == 3)
  error ('Not implemented yet')
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