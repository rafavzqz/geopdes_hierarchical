coeficientes = 1;

switch coeficientes
    case 1,
        ctes.gamma = 1/4;
        ctes.sigma = -7/4*pi;
    case 2, % genera coeficientes con un mayor salto
        ctes.gamma = 0.1;
        ctes.sigma = -14.92256510455152;
end

ctes.rho = pi/4;
ctes.R =-tan((pi/2-ctes.sigma)*ctes.gamma)/tan(ctes.rho*ctes.gamma);
ctes.a1 = ctes.R;
ctes.a2 = 1.0;
ctes.delta = 1/max(ctes.a1,ctes.a2);

% diffusion coefficient (a) of the equation
problem_data.diff = @(x) ctes.a1*(x(1)*x(2)>0)+ctes.a2*(x(1)*x(2)<0);

f = @(x, y) zeros (size (x));
% Dirichlet data, function g_D
h = @(x,y,ind) singular_solution(x,y,ctes.gamma,ctes.rho,ctes.sigma);


graduex = @(x) grad_singular_solution(x,ctes.gamma,ctes.rho,ctes.sigma);
uex = @(x,y) singular_solution(x,y,ctes.gamma,ctes.rho,ctes.sigma);


drchlt_sides = [1 2 3 4];
dim = 2;