function Animate(p)
files = dir([p,'/*.tif']);
for i = 1:length(files)
    f=[p, '/', files(i).name];
    im=imread(f);
    Nim=double(im)/255;
    FindRect(Nim);
    %drawnow;
    pause(0.5);
    
end
end

