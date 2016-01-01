%% Introduction:
% The code simmulates stretching a dsDNA molecule in optical tweezers.
% The goal is to fit WLC model to the simulated data (not in this code) to
% study the effect of averaging the data and sampling rate on the fitting
% parameters and chi squared value.

%TO DO::
% increase simulation rate
% look at the histogram of points at each extension (gaussian?)

clc
clear
cla
close all
%% parameters and values
% unit are in SI

nmperbp = 0.336; %nm/basepair of DNA
bp = 1000;
L = nmperbp*bp*1e-9; %DNA contour length (m)
P = 36.24e-9; %DNA persistence length (m)
fs = 5e5; %simulation frequency Hz

kb =1.3806488e-23;
temp = 298; %Kelvin
kT = kb*temp;
f_max = 11e-12; %N
dts = 1.0/fs; %simulation time in seconds
k = 100e-6; %trap stiffness N/m
m =5e4; %pipette dwelling counts
mm = 5e3; %save data every mm counts
bin_data = m/mm; %number of data saved in each pipette dwell

% calculate constants:
G=6*pi*1e-3*2.89*1e-6/2.0; %all units SI
fc=k/(2.0*pi*G);
D=kb*temp/G;
D1=(-1.0)*2.0*pi*fc; %position is included in the main formula
D2=2.0*D;

f_step = [];
ext_step = [];
x_pipette = [];
v_pipette = 300e-9; %m/sec pipette velocity

force = [];
force_std = [];
force_sem = [];
extension = [];
extension_std = [];
extension_sem = [];

n=1;

%% Stretching in OT

x_pipette(1) = 10.0e-9;
x_trap = 0.0;
ext = x_pipette(1)-x_trap;
WLC = ((kT/P)*(1/(2*(1 - (ext/L)))^2 + ...
    (ext/L) - (1/4)));
f= k*x_trap;

while f <f_max
    %next pipette step
    n=n+1;
    x_trap_old = x_trap;
    x_pipette(n) = x_pipette(n-1)+(v_pipette*dts);
    
    if mod(n,10000)==0
        u=0;
        
        %TODO:: find best values for m and sampling rate for uncorrelated
        %fluctuations
        for i = 1:m
            x_trap = x_trap_old*(1+(D1*dts))+ ...
                (WLC*dts)/G+ ...
                (sqrt(D2*dts)*normrnd(0.0, 1.0)); %langevin equation
            ext = x_pipette(n)-x_trap;
            WLC = ((kT/P)*(1/(2*(1 - (ext/L)))^2 + ...
                (ext/L) - (1/4)));
            f = k*x_trap;
            x_trap_old = x_trap;
            
            if mod(i,2500)==0 %ensures uncorrelated fluctuations
                u=u+1;
                ext_step(u) = ext;
                f_step(u) = f;
            end
            
        end
        
        %averaging extension and force at each pipette step
        
        force(end+1) = mean(f_step);
        force_std(end+1) = std(f_step);
        force_sem(end+1) = force_std(end)/sqrt(length(f_step));
        
        extension(end+1) = mean(ext_step);
        extension_std(end+1) = std(ext_step);
        extension_sem(end+1) = extension_std(end)/ ...
            sqrt(length(ext_step));
        
    else
        
        x_trap = x_trap_old*(1+(D1*dts))+ ...
            (WLC*dts)/G+ ...
            (sqrt(D2*dts)*normrnd(0.0, 1.0)); %langevin equ.
        ext = x_pipette(n)-x_trap;
        WLC = ((kT/P)*(1/(2*(1 - (ext/L)))^2 + ...
            (ext/L) - (1/4)));
        f = k*x_trap;
    end
    
    %display last results on the screen to see the progress
    if mod(n, 50000)==0
        sprintf ('%f' , force(end)*1e12)
        sprintf ('%f' , extension(end)*1e9)
    end
    
    
end

% % windowing data
% ext_sub = extension(1:end-mod(length(extension),w));
% ext_avg = mean(reshape(ext_sub,w,[]));
% ext_std = std(reshape(ext_sub,w,[]));
% %
% f_sub = force(1:end-mod(length(force),w));
% f_avg = mean(reshape(f_sub,w,[]));
% f_std = std(reshape(f_sub,w,[]));
% %
% Nf_std = f_std/sqrt(w);
%
% output_avg = [ext_avg'*1e9 f_avg'*1e12 f_std'*1e12 Nf_std'*1e12];

%% saving the results
output_avg = [extension'*1e9 force'*1e12 force_sem'*1e12];
%filename includes DNA info, sampling rate, number of data points in each 
%pipette dwell and dwell counts
filename = sprintf('DNA%dkbp_simRate%dkHz_%d_dwellcount%d.txt', ...
    bp/1000, fs/1000,bin_data,m); 
dlmwrite(filename, output_avg);

%% Results

f1=figure(1);
hold on; grid on;
set(gca, 'FontSize',16);
plot(extension*1e9,force*1e12);
name = sprintf('DNA: %d kbp,  simulation rate: %d kHz', bp/1000, fs/1000);
title(name);
xlabel('average extension(nm)');
ylabel('average force(pN)');

% f2=figure(2); hold on;
% grid on;
% plot(extension,force_std);
% name = sprintf('DNA: %d kbp,  simulation rate: %d kHz', bp/1000, fs/1000);
% title(name);
% xlabel('average extension(m)');
% ylabel('force standard deviation(N)');
%
% f3=figure(3); hold on;
% grid on;
% plot(extension,extension_std);
% name = sprintf('DNA: %d kbp,  simulation rate: %d kHz', bp/1000, fs/1000);
% title(name);
% xlabel('average extension(m)');
% ylabel('extension standard deviation(m)');
