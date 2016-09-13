clc
clear
cla

figure(1); hold on; box on;

maxN = 75;
cmap = jet(maxN+1);


for non0 = 0:1:maxN
    
    [time,Focal]=ode45(@nON, [0 9], non0);
    
    if(mod (non0,5)==0)
        
        plot (time,Focal, '-', 'Color', cmap(non0+1, :),'LineWidth', 1.5);
        
        drawnow;
    end
end

set(gca, 'FontSize', 16);
xlabel('time (sec)');
ylabel('number of focal adhesions');
title ('0 < non < 75  &  ntot=75');