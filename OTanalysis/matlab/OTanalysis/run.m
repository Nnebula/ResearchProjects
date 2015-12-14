close all
clc
clear

TrapFile = '212903trap_2.txt';
PipetteFile = '212903pipette_2.txt';
PSDfile = '2011-02-15-16-58-16_PSD.txt';

[ yTrap_um yPipette_um NyPSD yPiezo yForce_Image yExtension_Image ] = getData(TrapFile, PipetteFile, PSDfile );

% fZero_Image = findZero(yForce_Image)
% yForce_Image0 = yForce_Image - fZero_Image;
% plot(yExtension_Image, yForce_Image0);

cycleIndex = getCycle(yTrap_um);
yTrap_um_sub = yTrap_um(cycleIndex(1):cycleIndex(2));
NyPSD_sub = NyPSD(cycleIndex(1):cycleIndex(2));

NyPSD_um = MicronPerVolt(yTrap_um_sub, NyPSD_sub);

yExtension_PSD = yPiezo - NyPSD_um;
yForce_PSD = NyPSD_um * 113.2;

plot(yExtension_PSD,yForce_PSD)