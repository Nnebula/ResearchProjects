clear
d=load('SA_f0_20NSEG_ds005_sweep1000.csv');

figure(1)

for N = 1:length(d)
    % plot chain
    plot3(d(N, 1:3:end),d(N, 2:3:end), d(N, 3:3:end), ...
                '-b','LineWidth', 1,...
                'MarkerEdgeColor','k', ...
                'MarkerFaceColor','r', ...
                'MarkerSize', 5)
    hold all
    % plot start and end points
    plot3(d(N, 1),d(N, 2), d(N, 3), ...
                'o','LineWidth', 1,...
                'MarkerEdgeColor','k', ...
                'MarkerFaceColor','r', ...
                'MarkerSize', 5)
    plot3(d(N, end-2),d(N, end-1), d(N, end), ...
                'o','LineWidth', 1,...
                'MarkerEdgeColor','k', ...
                'MarkerFaceColor','r', ...
                'MarkerSize', 5)
    hold off
    grid on
    axis([-5 5 -5 5 -5 5])
    
    % create animation and save to animated gif
    drawnow
    filename = 'sachain05.gif';
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if N == 1;
        imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif', 'WriteMode','append', 'DelayTime', 0);
    end
end

