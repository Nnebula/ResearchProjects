clc
clear

files=dir('*.txt');

f1=figure()
for i=1:length(files)
    filename=files(i).name;
    data=load(filename, '-ascii');
    
    dataname=strrep(filename,'_','-');
    
    set(f1, 'Name', filename)
    subplot(2,1,1);
    plot(data(:,1) , data(:,2), 'r','LineWidth',1.5)
    grid on
    legend(['CD' '-' dataname],'Location','SouthWest')
    xlabel ('wavelength(nm)')
    ylabel ('CD(mdeg)')
    
    subplot(2,1,2);
    plot(data(:,1) , data(:,3), 'r','LineWidth',1.5)
    grid on
    legend(['HT' '-' dataname] ,'Location','SouthWest')
    xlabel ('wavelength(nm)')
    ylabel ('High Tension (V)')
    
    print(f1,'-dpng',[dataname '.png'])
end

close all