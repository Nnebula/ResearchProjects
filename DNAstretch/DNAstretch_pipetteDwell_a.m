%% Introduction:
% The code simmulates stretching a dsDNA molecule in optical tweezers.
% The goal is to fit WLC model to the simulated data (not in this code) to
% study the effect of averaging the data and sampling rate on the fitting
% parameters and chi squared value.

clc
clear
close all
%% parameters and values
% unit are in SI

nmperbp = 0.336; %nm/basepair of DNA
bp = 1000;
L = nmperbp*bp*1e-9; %DNA contour length (m)
P = 36.24e-9; %DNA persistence length (m)
fs = 100000; %simulation frequency Hz

kb =1.3806488e-23;
temp = 298;
kT = kb*temp;
f_max = 12.0e-12; %N
dts = 1.0/fs; %simulation time in seconds
k = 100e-6; %trap stiffness N/m
m =10000; %pipette dwelling counts

% calculate constants:
G=6*pi*1e-3*2.89*1e-6/2.0; %all units SI
fc=k/(2.0*pi*G);
D=kb*temp/G;
D1=(-1.0)*2.0*pi*fc; %position is included in the main formula
D2=2.0*D;

f_step = zeros(1, m);
ext_step = zeros(1, m);
x_pipette = [];
v_pipette = 100e-9; %m/sec pipette velocity

force = [];
force_std = [];
force_Nstd = [];
extension = [];
extension_std = [];
extension_Nstd = [];

n=1;
w = 20;

%% Stretching in OT

x_pipette(1) = 10.0e-9;
x_trap = 0.0;
ext = x_pipette(1)-x_trap;
WLC = ((kT/P)*(1/(2*(1 - (ext/L)))^2 + ...
    (ext/L) - (1/4)));
f= k*x_trap;
j=0;

while f <f_max
    %next pipette step
    n=n+1;
    x_trap_old = x_trap;
    x_pipette(n) = x_pipette(n-1)+(v_pipette*dts);
    
    if mod(n,2000)==0
        
        for i = 1:m
            x_trap = x_trap_old*(1+(D1*dts))+ ...
                (WLC*dts)/G+ ...
                (sqrt(D2*dts)*normrnd(0.0, 1.0)); %langevin equ.
            ext_step(i) = x_pipette(n)-x_trap;
            WLC = ((kT/P)*(1/(2*(1 - (ext_step(i)/L)))^2 + ...
                (ext_step(i)/L) - (1/4)));
            f_step(i) = k*x_trap;
            x_trap_old = x_trap;
        end
        %averaging extension and force at each pipette step
        j =j+1;
        
        force(j) = mean(f_step);
        force_std(j) = std(f_step);
        force_Nstd(j) = force_std(j)/sqrt(m);
        
        extension(j) = mean(ext_step);
        extension_std(j) = std(ext_step);
        extension_Nstd(j) = extension_std(j)/sqrt(m);
        
    else
        
        x_trap = x_trap_old*(1+(D1*dts))+ ...
            (WLC*dts)/G+ ...
            (sqrt(D2*dts)*normrnd(0.0, 1.0)); %langevin equ.
        ext = x_pipette(n)-x_trap;
        WLC = ((kT/P)*(1/(2*(1 - (ext/L)))^2 + ...
            (ext/L) - (1/4)));
        f = k*x_trap;
    end
    
    if mod(n, 10000)==0
        sprintf ('%f' , force(j-1)*1e12)
        sprintf ('%f' , extension(j-1)*1e9)
    end
    
    if mod(n, 50000)==0
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
output_avg = [extension'*1e9 extension_std'*1e9 extension_Nstd'* ...
    1e9 force'*1e12 force_Nstd'*1e12 force_std'*1e12];
dlmwrite('DNA1kbp_10000.txt', output_avg);

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

