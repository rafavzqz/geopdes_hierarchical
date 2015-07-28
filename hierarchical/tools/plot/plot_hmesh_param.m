function plot_hmesh_param(hmsh, nfig)
%
% function plot_hmesh_param(hmsh, nfig)
%
% This function plots the parametric hierarchical mesh
%
% INPUT:    hmsh:
%           nfig: integer (Optional)
%
% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% ATENCION: Completar la descripcion de esta funcion y mejorar
% 
%

if hmsh.ndim == 1
    disp('Aun no esta implementado un grafico de la malla en 1d'),
    return,
end

if nargin == 1
    figure
else
    close(figure(nfig))
    figure(nfig)
end

switch hmsh.ndim
    case 2, axis([0 1 0 1])
    case 3, axis([0 1 0 1 0 1])
end
axis square
axis off
hold on

switch hmsh.ndim
    case 2,
        for i = 1:hmsh.nel
            lev = hmsh.globnum_active(i,1);
            mx = hmsh.mesh_of_level(lev).nel_dir(1);
            indx = hmsh.globnum_active(i,2);
            x = [(indx-1)/mx indx/mx];
            my = hmsh.mesh_of_level(lev).nel_dir(2);
            indy = hmsh.globnum_active(i,3);
            y = [(indy-1)/my indy/my];
            abcisas = [x(1) x(2) x(2) x(1) x(1)];
            ordenadas = [y(1) y(1) y(2) y(2) y(1)];
            plot(abcisas,ordenadas,'k');
        end
    case 3,
        for i = 1:hmsh.nel
            lev = hmsh.globnum_active(i,1);
            mx = hmsh.mesh_of_level(lev).nel_dir(1);
            indx = hmsh.globnum_active(i,2);
            x = [(indx-1)/mx indx/mx];
            my = hmsh.mesh_of_level(lev).nel_dir(2);
            indy = hmsh.globnum_active(i,3);
            y = [(indy-1)/my indy/my];
            mz = hmsh.mesh_of_level(lev).nel_dir(3);
            indz = hmsh.globnum_active(i,4);
            z = [(indz-1)/mz indz/mz];
            abcisas = [x; x(2) x(2); x(2) x(1); x(1) x(1); ...
                x; x(2) x(2); x(2) x(1); x(1) x(1); ...
                x(1) x(1); x(2) x(2); x(2) x(2); x(1) x(1)]';
            ordenadas = [y(1) y(1); y; y(2) y(2); y(2) y(1); ...
                y(1) y(1); y; y(2) y(2); y(2) y(1); ...
                y(1) y(1); y(1) y(1); y(2) y(2); y(2) y(2)]';
            cotas = [z(1)*ones(4,2); ...
                z(2)*ones(4,2); ...
                z; z; z; z]';
            plot3(abcisas,ordenadas,cotas,'k');
        end
end
hold off