close all
clear all
warning ('off','geopdes:nrbmultipatch')
% PHYSICAL DATA OF THE PROBLEM

% 3 patch surf hyperbolic smooth boundary
coefs_1(:,1,:)=[-0.5 -0.5 0.0; -0.375 -0.5 -0.125; -0.25 -0.5 -0.1875];
coefs_1(:,2,:)=[-0.5 0.0 0.5; -0.3 -0.15 0.15; -0.1 -0.3 -0.0625];
coefs_1(:,3,:)=[-0.5 0.5 0.0; -0.225 0.2 0.025; 0.05 -0.1 -0.0075];

coefs_2(:,1,:)=[0.5 -0.5 0.0; 0.5 -0.225 0.275; 0.5 0.05 0.2475];
coefs_2(:,2,:)=[0.125 -0.5 -0.375; 0.2 -0.2625 -0.0625; 0.275 -0.025 0.03];
coefs_2(:,3,:)=[-0.25 -0.5 -0.1875; -0.1 -0.3 -0.0625; 0.05 -0.1 -0.0075];

coefs_3(:,1,:)=[0.05 -0.1 -0.0075; 0.275 -0.025 0.03; 0.5 0.05 0.2475];
coefs_3(:,2,:)=[-0.225 0.2 0.025; 0.1375 0.2375 -0.1; 0.5 0.275 0.225];
coefs_3(:,3,:)=[-0.5 0.5 0.0; 0.0 0.5 -0.5; 0.5 0.5 0.0];

knots{1} = [0 0 0 1 1 1];
knots{2} = [0 0 0 1 1 1];

butterf1 = nrbmak(permute(coefs_1,[3,2,1]),knots);
butterf2 = nrbmak(permute(coefs_2,[3,2,1]),knots);
butterf3 = nrbmak(permute(coefs_3,[3,2,1]),knots);

nrb(1)=butterf1;
nrb(2)=butterf2;
nrb(3)=butterf3;
problem_data.geo_name = nrb;

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];
problem_data.weak_drchlt_sides = [];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));
% Source and boundary terms
C = 1;
p=-1;

% hyperboloid exponential (0, 0, 0)
cons=200;
problem_data.f = @(x, y, z) 16.*cons.*exp(1).^((-1).*cons.*(x.^2+y.^2+(x.^2+(-1).*y.^2).^2)).* ...
                            (1+4.*x.^2+4.*y.^2).^(-5).*(cons.^3.*(1+4.*x.^2+4.*y.^2).^3.*(4.* ...
                            x.^6+(-4).*x.^4.*((-1)+y.^2)+(y+2.*y.^3).^2+x.^2.*(1+8.*y.^2+(-4) ...
                            .*y.^4)).^2+8.*(64.*x.^8+x.^6.*(32+(-896).*y.^2)+2.*x.^4.*(7+(-48) ...
                            .*y.^2+832.*y.^4)+y.^2.*((-1)+14.*y.^2+32.*y.^4+64.*y.^6)+(-1).* ...
                            x.^2.*(1+60.*y.^2+96.*y.^4+896.*y.^6))+(-4).*cons.^2.*(1+4.*x.^2+ ...
                            4.*y.^2).^2.*(72.*x.^10+20.*x.^8.*(5+2.*y.^2)+x.^6.*(50+272.*y.^2+ ...
                            (-112).*y.^4)+(y+2.*y.^3).^2.*(1+7.*y.^2+18.*y.^4)+x.^4.*(11+142.* ...
                            y.^2+280.*y.^4+(-112).*y.^6)+x.^2.*(1+26.*y.^2+142.*y.^4+272.* ...
                            y.^6+40.*y.^8))+2.*cons.*(1+1568.*x.^10+19.*y.^2+174.*y.^4+796.* ...
                            y.^6+1752.*y.^8+1568.*y.^10+8.*x.^8.*(219+404.*y.^2)+4.*x.^6.*( ...
                            199+1064.*y.^2+2896.*y.^4)+2.*x.^4.*(87+786.*y.^2+3720.*y.^4+ ...
                            5792.*y.^6)+x.^2.*(19+244.*y.^2+1572.*y.^4+4256.*y.^6+3232.*y.^8)));
problem_data.h = @(x, y, z, ind) exp(-cons*(x.^2 + y.^2 + (x.^2 - y.^2).^2));

problem_data.uex     = @(x, y, z) exp(-cons*(x.^2 + y.^2 + (x.^2 - y.^2).^2));
problem_data.graduex = @(x, y, z) cat (1, ...
                       reshape ( (-2).*cons.*exp(1).^((-1).*cons.*(x.^2+y.^2+(x.^2+(-1).*y.^2).^2)).*x.*(1+4.*x.^2+4.*y.^2).^(-1).*(1+2.*x.^2+6.*y.^2), [1, size(x)]), ...
                       reshape ( (-2).*cons.*exp(1).^((-1).*cons.*(x.^2+y.^2+(x.^2+(-1).*y.^2).^2)).*y.*(1+6.*x.^2+2.*y.^2).*(1+4.*x.^2+4.*y.^2).^(-1), [1, size(x)]), ...
                       reshape ( (-4).*cons.*exp (1).^((-1).*cons.*(x.^2+y.^2+(x.^2+(-1).*y.^2).^2)).*(1+4.*x.^2+4.*y.^2).^(-1).*(x.^2+2.*x.^4+(-1).*y.^2+(-2).*y.^4), [1, size(x)])); 
                   
problem_data.lapuex  = @(x, y, z) 4.*cons.*exp(1).^((-1).*cons.*(x.^2+y.^2+(x.^2+(-1).*y.^2).^2)).*( ...
                                    1+4.*x.^2+4.*y.^2).^(-2).*((-1)+20.*cons.*x.^6+16.*cons.*x.^8+(( ...
                                    -8)+cons).*y.^2+4.*((-5)+2.*cons).*y.^4+20.*cons.*y.^6+16.*cons.* ...
                                    y.^8+x.^2.*((-8)+cons+(-24).*y.^2+16.*cons.*y.^2+44.*cons.*y.^4)+( ...
                                    -4).*x.^4.*(5+cons.*((-2)+(-11).*y.^2+8.*y.^4)));


degrees = [3];
for deg = degrees
for unif=0:0
  for trunc=1:1
        savename_mat=sprintf('3p_hyperb_bilap_mpC1_degree%d_unif%dtrunc%d.mat',deg,unif,trunc);
        savename_mesh=sprintf('3p_hyperb_bilap_mpC1_degree%d_unif%dtrunc%d_mesh.jpg',deg,unif,trunc);
        savename_solutions=sprintf('3p_hyperb_bilap_mpC1_degree%d_unif%dtrunc%d_sols.jpg',deg,unif,trunc);
% CHOICE OF THE DISCRETIZATION PARAMETERS (Coarse mesh)
clear method_data
method_data.degree      = deg*[1 1];        % Degree of the splines
method_data.regularity  = (deg-2)*[1 1];        % Regularity of the splines
method_data.nsub_coarse = [4 4];        % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];        % Number of subdivisions for each refinement
method_data.nquad      = 2*method_data.degree+5;       % Points for the Gaussian quadrature rule
method_data.Cpen = 10 * (min(method_data.degree) + 1);
method_data.space_type  = 'standard'; % 'simplified' (only children functions) or 'standard' (full basis)
method_data.truncated   = trunc;            % 0: False, 1: True

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag = 'elements';
adaptivity_data.C0_est = 1.0;
adaptivity_data.mark_param = .25;
adaptivity_data.mark_strategy = 'GERS'; % GERS -> hierarchical, GR -> global refinement (uniform)
adaptivity_data.max_level = 5;
adaptivity_data.max_ndof = 200000;
if unif==0
    adaptivity_data.num_max_iter = 6;
else
    adaptivity_data.num_max_iter = adaptivity_data.max_level-1;
end
adaptivity_data.max_nel = 200000;
adaptivity_data.tol = 1e-8;
adaptivity_data.adm_class=3;
if trunc==0
    adaptivity_data.adm_type = 'H-admissible';
else
    adaptivity_data.adm_type = 'T-admissible';
end

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

if unif==0
    [geometry, hmsh, hspace, u, solution_data] = adaptivity_bilaplace_mp_C1 (problem_data, method_data, adaptivity_data, plot_data);
elseif trunc==0
    [geometry, hmsh, hspace, u, solution_data] = uniform_bilaplace_mp_C1 (problem_data, method_data, adaptivity_data, plot_data);
end
if unif~=1 || trunc ~=1
    save(savename_mat,'geometry', 'hmsh', 'hspace', 'u', 'solution_data');
end

if unif~=1
% EXPORT VTK FILE

subplot(1,3,2)
npts = [51 51];
sp_plot_solution (u, hspace, geometry, hmsh, npts); shading interp
vtk_pts = {linspace(0,1,npts(1)), linspace(0,1,npts(2))};
C = hspace_subdivision_matrix (hspace);
% for iptc = 1:hmsh.npatch
%   F = reshape (geometry(iptc).map(vtk_pts), [3, npts]);
%   X = reshape (F(1,:), npts); Y = reshape (F(2,:), npts); Z = reshape (F(3,:), npts);
%   subplot(1,3,3)  
%   surf(X, Y, Z, problem_data.uex(X,Y,Z));
%   hold on; shading interp
%   [val] = sp_eval (hspace.space_of_level(end).Cpatch{iptc}*C{end}*u, hspace.space_of_level(end).sp_patch{iptc}, geometry(iptc), [npts]);
%   subplot(1,3,2)
%   surf(X, Y, Z, problem_data.uex(X,Y,Z)-val);
%   hold on; shading interp
% end

subplot(1,3,1)
npts = [51 51];
sp_plot_solution (u, hspace, geometry, hmsh, npts); shading interp
saveas(gcf,savename_solutions);

% EXPORT VTK FILE
npts = [51 51];
output_file=['3p_hyperboloid_exp_non_uniform_degree' num2str(deg)];
sp_to_vtk (u, hspace, geometry, npts, output_file, {'solution', 'gradient', 'laplacian'}, {'value', 'gradient', 'laplacian'})

%MESH PLOT
figure(2)
for ii = 1:numel(geometry); nrbplot (geometry(ii).nurbs, [100 100]); 
hold on; end; shading interp
hmsh_plot_cells (hmsh, 11, 11);
figure (2)
view(35,38)
set(gcf,'color','white')
axis equal
axis off
saveas(gcf,savename_mesh);

end
end
end
end
warning ('on','geopdes:nrbmultipatch')