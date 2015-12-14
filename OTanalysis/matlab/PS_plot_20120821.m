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
    loglog(data(:,1) , data(:,3), 'r','LineWidth',1.5)
    grid on
    legend(['PSD' '-' dataname],'Location','NorthEast')
    xlabel ('Frequency(Hz)')
    ylabel ('PSD(V^2/Hz)')
    title ('Sine PSD')

    print(f1,'-dpng',[dataname '.png'])
end

close all