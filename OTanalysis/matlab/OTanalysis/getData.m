function [ yTrap_um yPipette_um NyPSD yPiezo yForce_Image yExtension_Image ] = getData( TrapFile, PipetteFile, PSDfile )
% Loads trap and pipette image analysis data and DAQ_timing data.
    
PSDdata = load(PSDfile);

TrapData = importdata (TrapFile, '\t', 1);
TrapData = TrapData.data;

PipetteData = importdata (PipetteFile, '\t', 1);
PipetteData = PipetteData.data;

% Matches the start time of PSD and Image file
tStart = max(PSDdata(1,1), TrapData(1,1)); 
ind_tPSD = find(PSDdata(:,1)==tStart, 1);
ind_tTrap = find(TrapData(:,1)==tStart, 1);

% Filter out non matching times at the start of the files
PSDdata = PSDdata(ind_tPSD:end, : );
TrapData = TrapData(ind_tTrap:end, : );
PipetteData = PipetteData(ind_tTrap:end, : );

NyPSD = PSDdata(:,3)./PSDdata(:,4);
yPiezo = PSDdata(:,5)*-4.8;
yTrap_um = TrapData(:,3)*0.0445;
yPipette_um = PipetteData(:,3)*0.0445;

yForce_Image = yTrap_um *113.2;
yExtension_Image = yPipette_um - yTrap_um;

% figure(1);
% plot(NyPSD, yTrap_um);
% 
% figure(2);
% plot(yPipette_um, yTrap_um);

end