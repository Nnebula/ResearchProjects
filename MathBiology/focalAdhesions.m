clc
clear
cla

figure(1); hold on; box on;

maxN = 75;
cmap = jet(maxN+1);

nSteadyState1 = zeros(maxN,1);
nSteadyState2 = zeros(maxN , 1);

for nTOT = 1:1:maxN

    nON = nONGenerator(nTOT);
    [time,Focal]=ode45(nON, [0 35], nTOT);
        
    plot (time,Focal, '-', 'Color', cmap(nTOT+1, :), 'LineWidth', 1.5);
    
    drawnow;
    
    nSteadyState1(nTOT) = Focal(end);
    
    nONback = nONbackGenerator(nTOT);
    
    if(nTOT>25)
    [time2,Focal2]=ode45(nONback, [0 100],0.5*Focal(end));
    nSteadyState2(nTOT) = Focal2(end);
    end
end

set(gca, 'FontSize', 16);
xlabel('time (sec)');
ylabel('number of focal adhesions');
title ('0 < ntot < 76  &  non0 = ntot');

%%

figure(2); clf; hold on; box on;
plot(1:1:maxN,nSteadyState1, '-r', 'LineWidth', 1.5);
plot(1:1:maxN,nSteadyState2, '-b', 'LineWidth', 1.5);
set(gca, 'FontSize', 16);
xlabel('total number of focal adhesions (ntot)');
ylabel('steady state number of attached focal adhesions (non)');
title ('bifurcation');
