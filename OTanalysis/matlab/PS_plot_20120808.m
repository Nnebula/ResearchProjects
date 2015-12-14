clc
clear

files=dir('*.dat');

f1=figure();
for i=1:length(files)
    filename=files(i).name;
    data=load(filename, '-ascii');
    
    pos=strfind(filename, '_PowSpec');
    originaldataname=filename(1:pos-1);
    dataname=strrep(originaldataname,'_','-');
    
    set(f1, 'Name', filename)
    subplot(3,1,1);
    loglog(data(:,1) , data(:,2), 'r','LineWidth',1.5)
    grid on
    legend(['xPSD' '-' dataname],'Location','SouthWest')
    xlabel ('Frequency(Hz)')
    ylabel ('PSD(V^2/Hz)')
    
    subplot(3,1,2);
    loglog(data(:,1) , data(:,3), 'r','LineWidth',1.5)
    grid on
    legend(['yPSD' '-' dataname] ,'Location','SouthWest')
    xlabel ('Frequency(Hz)')
    ylabel ('PSD(V^2/Hz)')
    
    subplot(3,1,3)
    loglog(data(:,1) , data(:,2), 'r', data(:,1) , data(:,3), 'b')
    legend(['xPSD' '-' dataname],['yPSD' '-' dataname] ,'Location','SouthWest')
    xlabel ('Frequency(Hz)')
    ylabel ('PSD(V^2/Hz)')
    print(f1,'-dpng',[dataname '.png']);
end

close all