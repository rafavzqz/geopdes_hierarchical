function hspace = refine_hierarchical_space(hmsh, hspace, M, flag, new_cells, boundary)
%
% function hspace = refine_hierarchical_space(hmsh, hspace, M, flag, new_cells, boundary)
%
% This function updates the information in hspace. The variable hmsh has
% already been updated by refine_hierarchical_mesh
% ATENCION: Voy a limpiar y acomodar esta funcion
%
%
% Input:        M:  (1 x msh.nlevels cell array),
%                   where marked{lev} is a matrix whose rows are the tensor-product indices
%                   of marked functions or elements of level lev, for lev =
%                   1:hmsh.nlevels
%               flag: 'functions' or 'elements'
%
%
% Output:   hspace updated
%
% This function uses:    prepare_to_refine_space
%                        compute_functions_to_remove
%                        update_active_functions
%

if nargin == 4
    boundary = true;
end

% Computation of indices of functions of level lev that will become
% nonactive when removing the functions or elements in M{lev}
M = compute_functions_to_remove(hmsh, hspace, M, flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Enlarging space_of_level and Proj, if it is needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(hspace.space_of_level) < hmsh.nlevels
    % We fill hspace.space_of_level(hmsh.nlevels)
    
    tic 
    disp('Enlarging space_of_level:')
    
    msh_level = hmsh.mesh_of_level(hmsh.nlevels);
    [knots,aaa] = kntrefine (hspace.space_of_level(hmsh.nlevels-1).knots, hmsh.nsub-1, hspace.degree, hspace.degree-1);
    
    hspace.space_of_level(hmsh.nlevels) = sp_bspline (knots, hspace.degree, msh_level);
    
    coarse_space = hspace.space_of_level(hmsh.nlevels-1).constructor (msh_level);
    
tempo = toc;
fprintf(' %f seconds\n', tempo);
end

if size(hspace.Proj,1) < (hmsh.nlevels - 1)
    % We fill hspace.Proj{hmsh.nlevels-1,:}
    tic 
    disp('Enlarging Proj:')
    
    Proj = cell(1,hmsh.ndim);
    geo_1d = geo_load (nrbline ([0 0], [1 0]));
    msh = hmsh.mesh_of_level(hmsh.nlevels);
    for idim = 1:hmsh.ndim
        % Construct the 1D mesh
        mesh_1d = msh_cartesian ({msh.breaks{idim}}, msh.qn{idim}, msh.qw{idim}, geo_1d);
        mesh_1d = msh_precompute (mesh_1d);
        
        % Evaluate the coarse space in the fine mesh
        sp_coarse = coarse_space.sp_univ(idim);
        sp_fine = hspace.space_of_level(hmsh.nlevels).sp_univ(idim);
        
        Mass = op_u_v (sp_fine, sp_fine, mesh_1d, ones (mesh_1d.nqn, mesh_1d.nel));
        Gram = op_u_v (sp_coarse, sp_fine, mesh_1d, ones (mesh_1d.nqn, mesh_1d.nel)); % Gram matrix
        Pr = sparse (sp_fine.ndof, sp_coarse.ndof);
        for ii = 1:sp_coarse.ndof
            ind = find(Gram(:,ii));
            Pr(ind,ii) = Mass(ind,ind) \ Gram(ind, ii);
        end
        %     P = M \ G;
        % Provisorio:
        Pr(abs(Pr)<1e-5) = 0; % Mejorar esta linea
        %Proj{idim} = Pr;
        hspace.Proj{hmsh.nlevels-1,idim} = Pr; 
    end
    
    tempo = toc;
fprintf(' %f seconds\n', tempo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic 
disp('Preparing space:')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guardamos la informacion de las celdas activas en E (por comodidad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ne = cumsum([0; hmsh.nel_per_level(:)]);
E = cell(hmsh.nlevels,1);
% El siguiente loop seguramente se puede evitar usando mat2cell
for lev = 1:hmsh.nlevels
    ind_e = (Ne(lev)+1):Ne(lev+1); 
    E{lev} = hmsh.globnum_active(ind_e, 2:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of active functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndof_per_level = zeros(1, hspace.nlevels); 
degree = hspace.degree;
[A, W, R] = prepare_to_refine_space(hspace);
tempo = toc;
fprintf(' %f seconds\n', tempo);

tic
disp('Updating active dofs:')

aux_ind_fun = [];

if (boundary && hmsh.ndim > 1)
    new_elem = new_cells.interior;
else
    new_elem = new_cells;
end

nel_dir = cell(hmsh.nlevels,1);
ndof_dir = cell(hmsh.nlevels,1);
for i = 1: hmsh.nlevels
    nel_dir{i} = hmsh.mesh_of_level(i).nel_dir;
    ndof_dir{i} = hspace.space_of_level(i).ndof_dir;
end

% Deactivation of functions which have to be removed and activation of the
% new ones
[A, W, R] = update_active_functions(A, W, R, E, hmsh.removed, new_elem, M, degree, nel_dir, hspace.Proj, ndof_dir, hspace.type);

hspace.nlevels = numel(A);

for lev = 1:hspace.nlevels
ndof_per_level(lev) = size(A{lev},1);
    aux_ind_fun = vertcat(aux_ind_fun, lev*ones(ndof_per_level(lev),1));
end

hspace.ndof_per_level = ndof_per_level;
hspace.ndof = sum(hspace.ndof_per_level);
hspace.globnum_active = cell2mat(A);
hspace.globnum_active = horzcat(aux_ind_fun,hspace.globnum_active);
hspace.coeff = cell2mat(W);
hspace.removed = R;

hspace.active = cell(hspace.nlevels,1);
for lev = 1:hspace.nlevels
    % Mejorar lo siguiente
    switch hmsh.ndim
        case 1, hspace.active{lev} = A{lev};
        case 2, hspace.active{lev} = sub2ind(hspace.space_of_level(lev).ndof_dir,A{lev}(:,1),A{lev}(:,2));
        case 3, hspace.active{lev} = sub2ind(hspace.space_of_level(lev).ndof_dir,A{lev}(:,1),A{lev}(:,2), A{lev}(:,3));
    end
    % hspace.active{lev} = sort(hspace.active{lev});
end
tempo = toc;
fprintf(' %f seconds\n', tempo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Updating C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing C:')
aux_active = hspace.active;
dif = hmsh.nlevels - hspace.nlevels;
if dif
    aux_active = vertcat(hspace.active,cell(dif,1));
end
ndof = zeros(1,hmsh.nlevels);
for l = 1:hmsh.nlevels
     ndof(l) = hspace.space_of_level(l).ndof;
end
hspace.C = compute_matrices_for_changing_basis(hmsh.nlevels, aux_active, ndof, hspace.Proj);
tempo = toc;
fprintf(' %f seconds\n', tempo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of sp_lev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Computing sp_lev:')
hspace.sp_lev = cell(hmsh.nlevels,1);
% Unefficient way. This is ok for the first working version
for ilev = 1 : hmsh.nlevels
   hspace.sp_lev{ilev} = sp_evaluate_element_list (hspace.space_of_level(ilev), hmsh.msh_lev{ilev}, 'gradient', true,'hessian', true);
end
tempo = toc;
fprintf(' %f seconds\n', tempo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (boundary && hmsh.ndim > 1)
    for iside = 1:2*hmsh.ndim
        %%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
        ind = setdiff (1:hmsh.ndim, ceil(iside/2));
        i = setdiff (1:hmsh.ndim, ind);
        if mod(iside,2) == 1
            boundary_ind = ones(hspace.nlevels,1);
        else
            boundary_ind = zeros(hspace.nlevels,1);
            for lev = 1:hspace.nlevels
                boundary_ind(lev) = hspace.space_of_level(lev).ndof_dir(i);
            end
        end
        M_boundary = cell(size(M));
        for lev = 1:numel(M)
            M_boundary{lev} = M{lev}(M{lev}(:,i) == boundary_ind(lev),ind);
        end
        hspace.boundary(iside) = refine_hierarchical_space(hmsh.boundary(iside), hspace.boundary(iside), ...
            M_boundary, 'functions', new_cells.boundary{iside}, false);
        % Now, we fill hspace.boundary(iside).dofs
        globnum_active_boundary = [hspace.boundary(iside).globnum_active(:,1:i) boundary_ind(hspace.boundary(iside).globnum_active(:,1)) ...
            hspace.boundary(iside).globnum_active(:,(i+1):end)];          
        [unos, hspace.boundary(iside).dofs] = ismember(globnum_active_boundary,hspace.globnum_active,'rows');
        if any(unos~=1)
            disp('Warning: Error when computing hspace.boundary().dofs')
            pause
        end
    end
else
    hspace.boundary = [];
end
