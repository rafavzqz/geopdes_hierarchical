function plot_msh(hmsh, iter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: This function works only in the parametric domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if hmsh.ndim ~= 2
%     disp('Hasta ahora, solo esta implementado el grafico de mallas en 2d')
%     return,
% end

%close(figure(1))

figure
switch hmsh.ndim
    case 2, axis([0 1 0 1])
    case 3, axis([0 1 0 1 0 1])
end
axis square
axis off
if nargin == 2
    title(sprintf('iter = %d  - %d active elements',iter, hmsh.nel),'FontSize',22)
else
    % title(sprintf('%d active elements', hmsh.nel))
end
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