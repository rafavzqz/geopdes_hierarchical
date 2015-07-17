% graficar singular solution

gamma = 1/4;
rho = pi/4;
sigma = -7/4*pi;

%gamma = .1; rho = pi/4; sigma= -14.92256510455152;

x=linspace(-1,1);
[X,Y] = meshgrid(x,x);

xx=X(:);
yy = Y(:);

zz = singular_solution([xx yy],gamma,rho,sigma);
Z = reshape(zz,100,100);
figure(1)
mesh(X,Y,Z)

grad_zz = grad_singular_solution([xx yy],gamma,rho,sigma);
grad_Z1 = reshape(grad_zz(:,1),100,100);
grad_Z2 = reshape(grad_zz(:,2),100,100);

figure(2)
mesh(X,Y,grad_Z1)
axis([-1 1 -1 1 -.3 0])

figure(3)
mesh(X,Y,grad_Z2)
axis([-1 1 -1 1 -.3 0])