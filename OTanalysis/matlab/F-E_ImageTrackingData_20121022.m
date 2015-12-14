clc
clear
cla

% Load the data from image txt data files
flist = dir('*pipette.txt');
pipette = importdata(flist.name,'\t',2);
xp = pipette.data(:,2);
yp = pipette.data(:,3);

flist = dir('*trap2.txt');
trap = importdata(flist.name,'\t',2);
xt = trap.data(:,2);
yt = trap.data(:,3);

% Pixel to micron conversion
xp = xp*0.0444;
yp = yp*0.0444;

xt = xt*0.0444;
yt = yt*0.0444;

% Calculate extension and force
yext = yp-yt;
yforce = yt*100;

% Plot force-extension from image data
h = figure;
plot(yext , yforce, 'r','LineWidth',1.5)
grid on
xlabel ('extension(um)')
ylabel ('force(pN)')

print(h,'-dpng','Figure1.png')








