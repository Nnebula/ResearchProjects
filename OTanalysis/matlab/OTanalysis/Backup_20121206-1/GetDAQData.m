function [ force exty ] = getDAQData( filename )
%getDAQData Reads a *DAQ.dat file and returns Force and Ext.y 
    data=load(filename);
    a=-5.1*data(:,2)./data(:,4);
    b=-5.1*data(:,3)./data(:,4);
    c=5.08*data(:,5);
    d=5.08*data(:,6);
    % extx = c + a;
    exty = -d - b;
    force = 100 * b;
end
