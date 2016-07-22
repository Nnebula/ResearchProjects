function fileread(p)

files = dir([p,'/*.tif']);

fig1 = figure(1);
widths = [];
heights = [];
objpos = [];

for i = 1:length(files)
    f=files(i).name;
    filename = strrep(f,'.tif', '');
    im=imread([p, '/', f]);
    Nim=double(im)/255;
    [w h found]=FindRect(Nim);
    
    if(found)
        heights(end+1)=h;
        widths(end+1)=w;
        objpos(end+1)=str2num(strrep(filename,'mV',''));
    end
    
    s=sum(Nim(400:480,:))/81;
    Nim2=Nim-ones(480,1)*s;
    
    sumx=sum(Nim2,1);
    sumy=sum(Nim2,2);
    
    set(fig1,'Name',filename);
    
    plot(sumy,'r','LineWidth',1.5 ); 
    hold on; plot( sumx,'b','LineWidth',1.5);
    hold off;
    
    set(gca,'FontSize', 18);
    
    xlabel ('position(pixel)')
    ylabel ('Intensity sum' , 'FontSize',18)
    legend('sum Y','sum X');
    
    print(fig1,'-dpng',[p, '/data/',filename '.png']);
    
end

f2=figure(2);
plot(objpos,heights, 'r', 'LineWidth',1.5 );
hold on;
plot(objpos, widths, 'b','LineWidth',1.5);
legend('Y size','X size');
hold off;

print(f2,'-dpng',[p, '/data/sizes.png']);

end
