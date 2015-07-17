function [problem_data, exact_solution] = load_problem_data(problem)

switch problem
    case 1,
        dospi = 2*pi;
        f = @(x,y) 2*dospi^2*sin(dospi*x).*sin(dospi*y);
        uex = @(x,y) sin(dospi*x).*sin(dospi*y);
        graduex = @(x,y) dospi*cat (1, ...
            reshape (cos(dospi*x).*sin(dospi*y), [1, size(x)]), ...
            reshape (sin(dospi*x).*cos(dospi*y), [1, size(x)]));
        dim = 2;
        drchlt_sides = 1:4;
    case 2,
        C = 40;
        normax2 = @(x,y) ((x-.5).^2+(y-.5).^2);
        uex = @(x,y) exp(-C*normax2(x,y));
        f = @(x,y) 4*C*(1-C*normax2(x,y)).*uex(x,y);
        graduex = @(x,y) -2*C*cat (1, ...
            reshape (uex(x,y).*(x-.5), [1, size(x)]), ...
            reshape (uex(x,y).*(y-.5), [1, size(x)]));
        dim = 2;
        drchlt_sides = 1:4;
    case 3,
        C = 100;
        normax2 = @(x,y) ((x-.5).^2+(y-.5).^2);
        uex = @(x,y) exp(-C*normax2(x,y));
        f = @(x,y) 4*C*(1-C*normax2(x,y)).*uex(x,y);
        graduex = @(x,y) -2*C*cat (1, ...
            reshape (uex(x,y).*(x-.5), [1, size(x)]), ...
            reshape (uex(x,y).*(y-.5), [1, size(x)]));
        
        dim = 2;
        drchlt_sides = 1:4;
    case 4;
        % Choose 1.5 < Cx, Cy < 2.5. Thus, the solution belongs H2\H3
        Cx = 1.6;
        Cy = 1.6;
        f = @(x,y) Cx*((Cx+1)*x.^(Cx-1)-(Cx-1)*x.^(Cx-2)).*y.^Cy.*(1-y)+...
            Cy*((Cy+1)*y.^(Cy-1)-(Cy-1)*y.^(Cy-2)).*x.^Cx.*(1-x);
        uex = @(x,y) x.^Cx.*(1-x).*y.^Cy.*(1-y);
        graduex = @(x,y) cat (1, ...
            reshape ((Cx*x.^(Cx-1).*(1-x)-x.^Cx).*y.^Cy.*(1-y), [1, size(x)]), ...
            reshape ((Cy*y.^(Cy-1).*(1-y)-y.^Cy).*x.^Cx.*(1-x), [1, size(x)]));
        dim = 2;
        drchlt_sides = 1:4;
    case 5;
        % Choose 1.5 < Cx, Cy < 2.5. Thus, the solution belongs H2\H3
        Cx = 2.1;
        Cy = 2.1;
        f = @(x,y) Cx*((Cx+1)*x.^(Cx-1)-(Cx-1)*x.^(Cx-2)).*y.^Cy.*(1-y)+...
            Cy*((Cy+1)*y.^(Cy-1)-(Cy-1)*y.^(Cy-2)).*x.^Cx.*(1-x);
        uex = @(x,y) x.^Cx.*(1-x).*y.^Cy.*(1-y);
        graduex = @(x,y) cat (1, ...
            reshape ((Cx*x.^(Cx-1).*(1-x)-x.^Cx).*y.^Cy.*(1-y), [1, size(x)]), ...
            reshape ((Cy*y.^(Cy-1).*(1-y)-y.^Cy).*x.^Cx.*(1-x), [1, size(x)]));
        dim = 2;
        drchlt_sides = 1:4;
    case 6;
        % Choose 1.5 < Cx, Cy < 2.5. Thus, the solution belongs H2\H3
        Cx = 1.6;
        Cy = 2.4;
        f = @(x,y) Cx*((Cx+1)*x.^(Cx-1)-(Cx-1)*x.^(Cx-2)).*y.^Cy.*(1-y)+...
            Cy*((Cy+1)*y.^(Cy-1)-(Cy-1)*y.^(Cy-2)).*x.^Cx.*(1-x);
        uex = @(x,y) x.^Cx.*(1-x).*y.^Cy.*(1-y);
        graduex = @(x,y) cat (1, ...
            reshape ((Cx*x.^(Cx-1).*(1-x)-x.^Cx).*y.^Cy.*(1-y), [1, size(x)]), ...
            reshape ((Cy*y.^(Cy-1).*(1-y)-y.^Cy).*x.^Cx.*(1-x), [1, size(x)]));
        dim = 2;
        drchlt_sides = 1:4;
    case 7,C = 100;
        normax2 = @(x,y,z) ((x-.5).^2+(y-.5).^2+(z-.5).^2);
        uex = @(x,y,z) exp(-C*normax2(x,y,z));
        f = @(x,y,z) 4*C*(6/4-C*normax2(x,y,z)).*uex(x,y,z);
        graduex = @(x,y,z) -2*C*cat (1, ...
            reshape (uex(x,y,z).*(x-.5), [1, size(x)]), ...
            reshape (uex(x,y,z).*(y-.5), [1, size(x)]), ...
            reshape (uex(x,y,z).*(z-.5), [1, size(x)]));
        dim = 3;
        drchlt_sides = 1:6;
    case 8,
        dospi = 2*pi;
        f = @(x,y,z) 3*dospi^2*sin(dospi*x).*sin(dospi*y).*sin(dospi*z);
        uex = @(x,y,z) sin(dospi*x).*sin(dospi*y).*sin(dospi*z);
        graduex = @(x,y,z) dospi*cat (1, ...
            reshape (cos(dospi*x).*sin(dospi*y).*sin(dospi*z), [1, size(x)]), ...
            reshape (sin(dospi*x).*cos(dospi*y).*sin(dospi*z), [1, size(x)]), ...
            reshape (sin(dospi*x).*sin(dospi*y).*cos(dospi*z), [1, size(x)]));
        dim = 3;
        drchlt_sides = 1:6;
    case 9, 
        C1 = 300; C2 = 300; C3=0; A1 = [7/12, 3/12]; A2 = [9/12,9/12]; A3 = [.75,.75];
        norma1x2 = @(x,y) ((x-A1(1)).^2+(y-A1(2)).^2);norma2x2 = @(x,y) ((x-A2(1)).^2+(y-A2(2)).^2);norma3x2 = @(x,y) ((x-A3(1)).^2+(y-A3(2)).^2);
        uex1 = @(x,y) exp(-C1*norma1x2(x,y));uex2 = @(x,y) exp(-C2*norma2x2(x,y));uex3 = @(x,y) exp(-C3*norma3x2(x,y));
        uex = @(x,y) uex1(x,y) + uex2(x,y) + uex3(x,y);
        f = @(x,y) 4*C1*(1-C1*norma1x2(x,y)).*uex1(x,y)+4*C2*(1-C2*norma2x2(x,y)).*uex2(x,y)+4*C3*(1-C3*norma3x2(x,y)).*uex3(x,y);
        graduex = @(x,y) -2*C1*cat (1, ...
            reshape (uex1(x,y).*(x-A1(1)), [1, size(x)]), ...
            reshape (uex1(x,y).*(y-A1(2)), [1, size(x)]))...
        -2*C2*cat (1, ...
            reshape (uex2(x,y).*(x-A2(1)), [1, size(x)]), ...
            reshape (uex2(x,y).*(y-A2(2)), [1, size(x)]))...
            -2*C3*cat (1, ...
            reshape (uex3(x,y).*(x-A3(1)), [1, size(x)]), ...
            reshape (uex3(x,y).*(y-A3(2)), [1, size(x)]));
        
        dim = 2;
        drchlt_sides = 1:6;
    case 10,
        dospi = 2*pi;
        f = @(x) dospi^2*sin(dospi*x);
        uex = @(x) sin(dospi*x);
        graduex = @(x) dospi*cos(dospi*x);
        dim = 1;
        drchlt_sides = 1:2;
    case 11,
        f = @(x, y) zeros (size (x));
        g = @test_square_g_nmnn;
        h = @(x, y, ind) exp (x) .* sin(y);
        
        uex     = @(x, y) exp (x) .* sin (y);
        graduex = @(x, y) cat (1, ...
            reshape (exp(x).*sin(y), [1, size(x)]), ...
            reshape (exp(x).*cos(y), [1, size(x)]));
        
        drchlt_sides = [1 2 3 4];
        dim = 2;
    case 12,
        chessboard_data
end

% dospi = 2*pi;
% switch dim
%     case 1, %f = @(x) sin(dospi*x);
%         f = @(x) x.*x;
%     case 2, f = @(x,y) sin(dospi*x).*sin(dospi*y);
%         %f = @(x,y) (1-x).^2.*(1-y).^2;
%     case 3, f = @(x,y,z) sin(dospi*x).*sin(dospi*y).*sin(dospi*z);
% end
h = @(x, y, ind) zeros(size(x));

switch dim
    case 1, geometry = geo_load (nrbline ([0 0], [1 0]));
    case 2, geometry = geo_load ('geo_square.txt'); %geometry=geo_load ('geo_ring.txt'); 
    case 3, geometry = geo_load ('geo_cube.txt');
end

nmnn_sides = setdiff(1:2*dim,drchlt_sides);

problem_data.geo_name = geometry;
problem_data.nmnn_sides = nmnn_sides;
problem_data.drchlt_sides = drchlt_sides;
%problem_data.c_diff = c_diff;
problem_data.f = f;
%problem_data.g = g;
problem_data.h = h;
problem_data.dim = dim;

exact_solution.uex = uex;
exact_solution.graduex = graduex;

