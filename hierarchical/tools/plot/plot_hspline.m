function plot_hspline(u, hmsh, hspace) 

% u degrees of freedom

[Z, pts] = hspline_eval(u, hmsh, hspace, 3);

npts = size(pts,2);
npts_idir = round(sqrt(npts));
figure(3)
for el = 1:hmsh.nel
    x = reshape(pts(1,:,el), npts_idir, npts_idir);
    y = reshape(pts(2,:,el), npts_idir, npts_idir);
    z = reshape(Z(:,el), npts_idir, npts_idir);
    surf(x,y,z);
    hold on
end
axis tight
xlabel('x')
ylabel('y')
hold off

% x = hmsh.geo_map(1,:,:);
% y = hmsh.geo_map(2,:,:);
% x = x(:); y = y(:);
% xxx = linspace(.01,.99,20);
% [X,Y] = meshgrid(xxx,xxx);
% ZZ = griddata(x,y,Z(:),X,Y);
% 
% surf(X,Y, ZZ);
% axis tight
% xlabel('x')
% ylabel('y')
