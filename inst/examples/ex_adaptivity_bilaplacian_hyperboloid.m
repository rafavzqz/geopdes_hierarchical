% PHYSICAL DATA OF THE PROBLEM
clear problem_data

% PHYSICAL DATA OF THE PROBLEM
problem_data.geo_name = 'geo_hyperboloid_ASG1.txt';

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4 5 6];
problem_data.weak_drchlt_sides = [];

% Physical parameters
problem_data.c_diff  = @(x, y, z) ones(size(x));

% Source and boundary terms
% Source and boundary terms
C = 1;
p = -1;
cons = 200;
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

% Exact solution
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

% DISCRETIZATION PARAMETERS
clear method_data
deg = 4;
method_data.degree      = [deg deg];     % Degree of the splines
method_data.regularity  = [deg-2 deg-2]; % Regularity of the splines
method_data.nsub_coarse = [4 4];         % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nsub_refine = [2 2];         % Number of subdivisions of the coarsest mesh, with respect to the mesh in geometry
method_data.nquad       = [deg+1 deg+1]; % Points for the Gaussian quadrature rule
method_data.space_type  = 'standard';
method_data.truncated   = 1;            % 0: False, 1: True
method_data.interface_regularity = 1;

% ADAPTIVITY PARAMETERS
clear adaptivity_data
adaptivity_data.flag          = 'elements';
adaptivity_data.C0_est        = 1.0;
adaptivity_data.mark_param    = .25;
adaptivity_data.mark_strategy = 'GERS';
adaptivity_data.max_level     = 10;
adaptivity_data.max_ndof      = 5000;
adaptivity_data.num_max_iter  = 8;
adaptivity_data.max_nel       = 5000;
adaptivity_data.tol           = 1e-5;

% GRAPHICS
plot_data.print_info = true;
plot_data.plot_hmesh = false;
plot_data.plot_discrete_sol = false;

[geometry, hmsh, hspace, u, solution_data] = adaptivity_bilaplace_mp_C1 (problem_data, method_data, adaptivity_data, plot_data);
npts = [21 21];
sp_plot_solution (u, hspace, geometry, npts); shading interp;
