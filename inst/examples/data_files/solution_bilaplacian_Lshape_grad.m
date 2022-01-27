

function graduex = solution_bilaplacian_Lshape_grad (x, y)
graduex = zeros(2, size(x,1), size(x,2));
  [th, r] = cart2pol (x, y);
  th = (th < 0).*(2*acos(-1) + th) + (th >= 0) .* th;
%f = @(z) sin(3*pi/2*z).^2 - z.^2 * sin(3*pi/2)^2;
z = 0.544483736782464; %z = fsolve(f, 0.5) The one on the left is more precise, computed with bisection method
% z = 0.908529189846099;

C1 = 1/(z-1)*sin(3*pi/2*(z-1)) - 1/(z+1)*sin(3*pi/2*(z+1));
C2 = cos(3*pi/2*(z-1)) - cos(3*pi/2*(z+1));
F1 = cos((z-1)*th) - cos((z+1)*th);
F2 = 1/(z-1)*sin((z-1)*th) - 1/(z+1)*sin((z+1)*th);


psi = C1*F1 - C2*F2;  %C1 and C2 are constants
Dpsi = C1*(-(z-1)*sin((z-1)*th)+(z+1)*sin((z+1)*th)) - C2*(cos((z-1)*th)-cos((z+1)*th));

graduex(1,:,:) = (z+1)*r.^(z).* psi.*cos(th) - r.^z.*Dpsi.*sin(th); 
graduex(2,:,:) = (z+1)*r.^(z).* psi.*sin(th) + r.^z.*Dpsi.*cos(th);
%   theta = (theta < 0).*(2*acos(-1) + theta) + (theta >= 0) .* theta;
%   f = rad.^(k*2/3) .* sin((k*2/3)*theta);

end

