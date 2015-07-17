function z = singular_solution(x,y,gamma,rho,sigma)


arg02pi = @(x,y)(mod(2*pi+atan2(y,x),2*pi));
normax2 = @(x,y)(x.^2+y.^2);

z = sqrt(normax2(x,y)).^gamma;

theta = arg02pi(x,y);

c = cos(gamma*[pi/2-sigma; rho; sigma; pi/2-rho]);

var = cos(gamma*[theta-pi/2+rho,theta-pi+sigma,theta-pi-rho,theta-3*pi/2-sigma]); 

var = var.*[0<=theta & theta< pi/2,pi/2<=theta & theta< pi, pi<=theta & theta< 3*pi/2, 3*pi/2<=theta & theta< 2*pi];

z = z.*(var*c);
