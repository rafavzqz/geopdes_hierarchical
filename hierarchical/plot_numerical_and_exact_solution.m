function plot_numerical_and_exact_solution(u, hmsh, hspace, uex)
%
% function plot_numerical_and_exact_solution(u, hmsh, hspace, uex)
%
% Funcion personal
% 
%

switch hmsh.ndim
    case 1, %x = hmsh.quad_nodes(:);
        [Z, x] = hspline_eval(u, hmsh, hspace, 4);
        figure(2)
        subplot (1,2,1)
        plot (x(:), Z(:),'*')
        title ('Numerical solution'), axis tight
        subplot (1,2,2)
        x= linspace(0,1);
        u_ex_values = uex(x);
        plot (x, u_ex_values)
        title ('Exact solution'), axis tight
        
    case 2,
        [Z, pts] = hspline_eval(u, hmsh, hspace, 3);        
        npts = size(pts,2);
        npts_idir = round(sqrt(npts));
        figure(2)
        subplot (1,2,1)
        for el = 1:hmsh.nel
            x = reshape(pts(1,:,el), npts_idir, npts_idir);
            y = reshape(pts(2,:,el), npts_idir, npts_idir);
            z = reshape(Z(:,el), npts_idir, npts_idir);
            surf(x,y,z);
            hold on
        end
        title ('Numerical solution'),
        axis tight
        xlabel('x')
        ylabel('y')
        hold off
        subplot (1,2,2)
        for el = 1:hmsh.nel
            x = reshape(pts(1,:,el), npts_idir, npts_idir);
            y = reshape(pts(2,:,el), npts_idir, npts_idir);
            z = uex(x,y);
            surf(x,y,z);
            hold on
        end
        title('Exact solution'),
        axis tight
        xlabel('x')
        ylabel('y')
        hold off
    case 3, disp('Todavia no esta implementado el grafico de la solucion discreta en 3d'),
end
