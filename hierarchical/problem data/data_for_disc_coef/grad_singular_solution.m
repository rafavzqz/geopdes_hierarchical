function y = grad_singular_solution(x,gamma,rho,sigma)


arg02pi = @(x)(mod(2*pi+atan2(x(:,2),x(:,1)),2*pi));
normax2 = @(x)(x(:,1).^2+x(:,2).^2);

y = gamma*sqrt(normax2(x)).^(gamma-2);

theta = arg02pi(x);

c = cos(gamma*[pi/2-sigma; rho; sigma; pi/2-rho]);

values = gamma*[theta-pi/2+rho,theta-pi+sigma,theta-pi-rho,theta-3*pi/2-sigma];

theta_matrix = [0<theta & theta< pi/2,pi/2<theta & theta< pi, pi<theta & theta< 3*pi/2, 3*pi/2<theta & theta< 2*pi];
 
mu = (cos(values).*theta_matrix)*c;

der_mu_sobre_gamma = -(sin(values).*theta_matrix)*c;

y = y*[1 1].*(mu*[1 1].*x+der_mu_sobre_gamma*[1 1].*[-x(:,2) x(:,1)]);

y(x(:,1)==0 | x(:,2)==0,:) = NaN;

